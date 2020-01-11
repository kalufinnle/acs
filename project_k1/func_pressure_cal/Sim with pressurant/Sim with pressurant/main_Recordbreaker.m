%--------------------------------------------------------------------------
%               UTAT ROCKETRY: DEFIANCE RECORDBREAKER
%                  HYBRID ROCKET ENGINE SIMULATION
%--------------------------------------------------------------------------

%Main script for simulating full rocket performance. Uses engine simulation
%and 2D kinematics simulation to determine altitude of rocket. Employs
%Euler's Method to step variables through time. Uses pre-combustion values 
%in kinematics model to provide conservative estimate. Ignores ignition 
%transient and throat erosion. Plots monitored variables against time.

%Original scripts: Thomas S.H. Leung
%Revisions: Mitchell L. Passarelli, Jacob Weber
%Version Date: Saturday, December 14, 2019

%%%%%GENRAL NOTES:
%%%%%   In Thomas' file "Project14.m" I use full line comments to track
%%%%%   lines I have either accounted for or do not need in this file.

%%%%%SHORT-TERM TO-DO:
%%%%%   1.) See if I can use analytical expression for port diameter to
%%%%%   find burn time (by solving d_i - d(t) = 0) to get length of
%%%%%   matrices of monitored variables. Also see if I can couple it with
%%%%%   ox_tank model and pressure calculator to find the minumum of the
%%%%%   three burnout conditions.
%%%%%   2.) Implement in house mass flow model into the ox tank sim
%%%%%   function.
%%%%%   3.) Nitrogen blowdown from pressurant tank.

%%%%%LONG-TERM TO-DO (i.e., after finished):
%%%%%   1.) Look for a way to vectorize the loop (start with engine calcs).
%%%%%       Looks like only the ox tank model and fuel regression need to
%%%%%       be looped.
%%%%%   2.) Implement a check for choked flow and execute different
%%%%%       functions for choked & unchoked flow.
%%%%%   3.) Implement CEA and coolprop wrappers for combustion and chemical
%%%%%       state properties.
%%%%%   4.) Implement vapour burn.
%%%%%   5.) Implement cross-wind wind effects.


clear; clc; close all;

%% Import Inputs From Excel & Initialize Variables %%
%importing values into imported
imported = import_params('Preliminary Engine Design.xlsm', 'Import Data_Recordbreaker', 1, 28);

%General parameters:
its = imported(1);          %Est'd # iterations.
dt = imported(2);           %Time step size (s).
theta = imported(3);        %Launch angle, from vertical (deg).
rail_len = 10;              %Launch rail length (m).
comb_eff = imported(23);    %Combustive efficiency.
nozz_eff = imported(24);    %Nozzle efficiency.
Pa = imported(4);           %Ambient pressure at launch altitude (Pa).
T_a = imported(5);          %Ambient temperature at launch altitude (K).
m = imported(6);            %Rocket dry mass (kg).
rho_fuel = imported(8);     %Mass density of fuel (kg/m^3).
rocket_length = imported(25);  %Total length of rocket (m).
launch_alt = imported(26); %used to calculate atmospheric density and ...

%CEA interpolation data:
num_rows = 525;     %# rows to import in CEA data file (= last row # - 1).
data = zeros(num_rows, 4);  %Cols. = OF, T, MW, gamma respectively.
[~, data(:, 1), data(:, 2), data(:, 3), ~, data(:, 4)] = ...
    import_cea_data('CEA_data.csv', 2, num_rows + 1);   %Get CEA data pts.
CEA = zeros(size(unique(data(:, 1)), 1), 4);    %Actual interp. pts matrix.
CEA(:, 1) = unique(data(:, 1));                 %Get unique OF values.

%Loop to clean up CEA data (average for each OF ratio):
for i = 1:size(CEA, 1)
    CEA(i, :) = mean(data(data(:, 1) == CEA(i, 1), :));
end
%NOTE: Don't need this loop if you use interp2 with OF and Pcc.

%Initialize variables and constants:
i = 1;              %Loop index reset.
dP_tol = 10;        % Presure drop tolerance per iteration  [Pa]
g0 = 9.81;          %Acceleration due to gravity at sea-level (m/s^2).
R_u = 8314;         %Universal gas constant (J/kmol*K).
MW_N2O = 44.013;    %Molecular weight/mass of N2O (kg/kmol).
MW_N2 = 28.0134;    % Molar mass of nitrogen [kg/kmol]
v_D = 0;            %Drag [velocity] loss (m/s).
f_D = 0.03;         %Darcy friction factor of flow sys.
L_fs = 0.00508;     %flow sys length (m).
d_fs = 2*0.027781;      %flow sys hydraulic diameter (m).
A_inj = imported(21); %injector area (m^2)
Cd = imported(22);  % discharge coefficient
%initialize variable arrays
s = zeros(2, its);  %Position vector (m).
v = s;              %Velocity vector (m/s).
a = s;              %Acceleration vector (m/s^2).
F = zeros(its, 1);  %Thrust (N). 
Total_Impulse = 0;  %Total Impulse (kNs).
Isp = F;            %Specific impulse (s).
OF = F;             %OF ratio.
Pcc = F;            %Combustion chamber pressure (Pa).
Pcc(1) = Pa(i);
Tcc = F;            %Combustion chamber temperature (K).
m_dot_ox = F;       %Ox mass flow rate (kg/s).
m_dot_f = F;        %Fuel mass flow rate (kg/s).
P_ox = F;           %Ox tank pressure array (Pa).
n_ox = zeros(its, 2);      % Oxidizer amount
C_D = F;            %Drag coefficient array
F_D = F;            %Drag force array (N).
m_f = F;            %Fuel mass array (kg).
M = F;              %Mach number.
m_f(1) = imported(7);       %Initial fuel mass (kg).
diam = imported(20);        %Diameter of rocket (m).
dtheta = 0;         %change in the angle of rocket with respect to vertical 
                    %this is zero since we are not accoutning for crosswind yet    
I_total = 0;        %Monitor total impulse. 
v_norm = F;         %Velocity magnitude.
a_norm = F;         %Acceleration magnitude.
compressive_force = F;  %Thrust added to drag.
G_ox = F;           %Oxidizer Flux
rho = F;            %Air Density
P_a = F;            %Atmospheric Pressure
P_a(1) = Pa;
P_dyn = F;          %Dynamic Pressure

% Pressurant tank parameters
n_pr_store = F;                 % Nitrogen amount 
n_pr = 0;                       % Initial pressurant amount (iterator)
n_pr_store(1, 1) = 0;           % Initial pressurant amount
n_pr_inc = 0.00000005;           % Pressurant increment (kmol)
%initial parameters
P_N2 = 31026407.82;         %Initial pressurant tank pressure (Pa) 4500psi.
R_N2 = R_u/MW_N2;
T_N2 = T_a;
vol_pr = 0.00147484;     %90 cu in pressurant tank (m^3). 
m_pr = vol_pr/(R_N2*T_N2/P_N2);
n_pr_avail = m_pr/MW_N2;    %amount of pressurant available
P_pr = F;                   % pressure in pressurant tank
P_pr(1) = P_N2;              

%Ox tank parameters:
pr_flg = 1; %-> Nitrogen
m_ox = F;
m_ox(1) = imported(9);      %Loaded liquid ox mass in tank (kg).
T_ox = F;
T_ox(1) = T_a;                 %Initial ox tank temperature (K).
V_tank = imported(10);      %Internal volume of ox tank (m^3).
m_tank = imported(11);      %Dry mass of ox tank (kg).
P_ox(1) = imported(28);     %Initial ox tank pressure (Pa).
Q = [2.781; 0.27244; 309.57; 0.2882];
        %Coefficients for molar specific vol. of liq. N2O (m^3/kmol).
Vhat_l = Q(2)^(1 + (1 - T_ox(1)/Q(3))^Q(4))/Q(1);
        %Initial molar specific vol. of liq. N2O (m^3/kmol).
n_ox(1,1) = m_ox(1)/MW_N2O;       %# moles of liquid ox; n(2) is gaseous moles.
n_tot = n_ox(1,1) + P_ox(1)*(V_tank - Vhat_l*n_ox(1,1))/(R_u*T_ox(1));
n_ox(1,2) = P_ox(1)*(V_tank - n_tot*Vhat_l)/(R_u*T_ox(1) - P_ox(1)*Vhat_l);
        %NOTE: Assuming P_sat = P_ox(1).

%Fuel core parameters:
a0 = imported(12);          %Regression rate coefficient (m/s^2).
n_reg = imported(13);       %Regression rate exponent.
num_ports = imported(14);   %# combustion ports in fuel grain.
L = imported(15);           %Fuel core length (m).
w = imported(16);           %Initial fuel web thickness (m).
r_fo = imported(17);        %Inner combustion chamber radius (m).
A_cc = pi*r_fo^2;           %Combustion chamber cross-section area (m^2).

%Nozzle parameters:
A_t = imported(18);         %[Initial] throat area (m^2).
A_e = imported(19);         %Exit area (m^2).
e = @(Ma_e, k) (2/(k+1)*(1 + (k-1)/2*Ma_e^2))^((k+1)/(2*k-2))/Ma_e;
        %Function handle for nozzle area expansion ratio.
% e = @(Ma_cc, k) A_cc/A_e*Ma_cc/(1 + (k-1)/2*Ma_cc^2)^((k+1)/(2*k-2));
%         %Function handle for chamber-exit area-Mach relation.
nps = imported(27); %pipe for rocket
addpath('Drag Model')
Drag_Model = Recordbreaker_AeroLab_Drag_Model(nps);


%% Determine whether flow sys or injector plate primarily constricts flow:
%%%%%NOTE: Don't worry about this one yet. Realistically, the injector
%%%%%plate should always constrict the flow.


%% Simulation Time Loop (ignition to apogee) %%
%Loop until rocket reaches apogee (ignore condition on first iteration):
while v(2, i) > 0 || i == 1
    is_burning = w > 0 && Pcc(i) <= P_ox(i) && m_ox(i) > 0;
    
    %If burn not finished, continue engine simulation:
    if is_burning
        burnout_alt = s(2, i);    %Update burnout altitude.
        burn_step = i;            %Step when engine finishes its burn
        %Simulate blow-down in ox tank:
        if i == 1   %Uses initial guesses for initial Pcc and m_dot_ox.
            [m_dot_ox(i), dn, dT_ox, P_ox(i + 1)] = ox_tank_sim(...
                n_ox(i, :), T_ox(i), V_tank, m_tank, 1.8, f_D, L_fs, d_fs, ...
                A_inj, Cd, 2.75e6, n_pr, pr_flg);
        else        %Uses prev. iteration data.
            [m_dot_ox(i), dn, dT_ox, P_ox(i + 1)] = ox_tank_sim(...
                n_ox(i, :), T_ox(i), V_tank, m_tank, m_dot_ox(i - 1), f_D, ...
                L_fs, d_fs, A_inj, Cd, Pcc(i - 1), n_pr, pr_flg);   
        end
       
        % Pressure drop between itterations
        dP = P_ox(i) - P_ox(i+1);
        %Add nitrogen to tank until pressure drop is minimized to within dp_tol
        while (dP>dP_tol)
            if (n_pr>=n_pr_avail)
                break
            end
            % Add nitrogen
            n_pr = n_pr + n_pr_inc;
            
            % update nitrogen tank pressure
            % modify the blowdown code?
            % going to have to have a check when the nitrogen tank pressure
            % drops below the nominal Pox
            
            % Recalculate blowdown
            [m_dot_ox(i), dn, dT_ox, P_ox(i + 1)] = ox_tank_sim(...
                    n_ox(i, :), T_ox(i), V_tank, m_tank, m_dot_ox(i - 1), f_D, ...
                    L_fs, d_fs, A_inj, Cd, Pcc(i - 1), n_pr, pr_flg); 

            % Calculate new pressure drop                                      
            dP = P_ox(i) - P_ox(i+1);          
        end
        
        
        %Determine mass flow rates and fuel regression:
        r_port = r_fo - w;                  %Get fuel port radius.
        G_ox(i) = m_dot_ox(i)/(pi*r_port^2);   %Calc. ox mass flux.
        reg_rate = a0*G_ox(i)^n_reg;           %Calc. fuel regression rate.
        SA_port = 2*pi*r_port*L;            %Calc. fuel core surface area.
        m_dot_f(i) = reg_rate*rho_fuel*SA_port*num_ports;
        m_dot_prop = m_dot_ox(i) + m_dot_f(i);  %Total mass flow rate.
        
        %Determine propellant thermochemical properties:
        OF(i) = m_dot_ox(i)/m_dot_f(i);
        Tcc(i) = interp1(CEA(:,1), CEA(:,2), OF(i));
        gamma = interp1(CEA(:,1), CEA(:,4), OF(i));
        R = R_u/interp1(CEA(:,1), CEA(:,3), OF(i));
        char_vel = comb_eff*sqrt(gamma*R*Tcc(i))/gamma * ...
            ((gamma + 1)/2)^((gamma + 1)/(2*gamma - 2));
        
        %Eval. engine performance including nozzle:
        if i > 1   %Calc. combustion chamber P (avg to damp out osc.).
            Pcc(i) = (m_dot_prop*char_vel/A_t + Pcc(i - 1))/2;
        else
            Pcc(i) = m_dot_prop*char_vel/A_t;
        end
        Ma_e = fzero(@(Ma) A_e/A_t - e(Ma, gamma), 2.86);   %Exit Mach.
        P_e = Pcc(i)/(1 + (gamma-1)/2*Ma_e^2)^(gamma/(gamma-1));
                %Flow is always choked for Pcc >= 200 psi.
        
%         %Eval. engine perfor mance when flow not choked:
%         rho = Pcc(i)/R/Tcc(i);              %Calc. chamber gas density.
%         v_cc = m_dot_prop/rho/A_cc;         %Get flow velocity.
%         Ma_cc = v_cc/sqrt(gamma*R*Tcc(i));  %Combustion chamber Mach #.
%         Ma_e = fzero(@(Ma_e) Ma_e/(1 + (gamma-1)/2*Ma_e)^...
%             ((gamma+1)/(2*gamma-2)) - e(Ma_cc, gamma), 0.5);    %Exit Mach.
        %%%%%Don't worry about above until sim is finished: can implement a
        %%%%%check for choked condition: look at AER310 midterm Q3 to see
        %%%%%how to check for choking.
        
        v_e = sqrt(2*gamma*R*Tcc(i)/(gamma - 1) * ...
            (1 - (P_e/Pcc(i))^((gamma - 1)/gamma)));    %Get exit velocity.
        
        if i == 1
            F(i) = nozz_eff*(m_dot_prop*v_e + (P_e - Pa)*A_e);     %Thrust.
        else
            F(i) = nozz_eff*(m_dot_prop*v_e + (P_e - P_a(i-1))*A_e);     %Thrust.
        end
       
        I_total = I_total + F(i)*dt;        %Tracking total impulse.
        Isp(i) = F(i)/m_dot_prop/g0;        %Calc. specific impulse.
        
        %Update engine parameters:
        %Make sure to store all parameters in (i+1)th entry, but pass only
        %the ith entries to the kinematics model.
        w = w - reg_rate*dt;                %Regress fuel core.
        m_f(i + 1) = m_f(i) - m_dot_f(i)*dt;    %Consume fuel mass.        
        m_ox(i + 1) = m_ox(i) + dn(1)*MW_N2O*dt; %Consume ox mass.
        n_ox(i+1, :) = n_ox(i, :) + dn*dt;                      %Consume moles of ox.
        T_ox(i + 1) = T_ox(i) + dT_ox*dt;          %Step ox tank temp.
        n_pr_store(i, 1) = n_pr;
        %%%%Include updates to ox tank parameters here too.%%%%
 
    %Change here for tracking
    else
        Pcc(i) = P_a(i-1);
        P_ox(i) = P_a(i-1);
        T_ox(i + 1) = T_ox(i);
        m_f(i + 1) = m_f(i);
        m_ox(i + 1) = m_ox(i)*(m_f(i) > 0);
    end

    %Move rocket forward in time using kinematics model:
    %%%%Call to kinematics model. Note that this needs to output the
    %%%%ambient pressure at the current altitude into P_a.
    %%%%Pull apart Thomas' KinematicsP14 function into atomic blocks and
    %%%%code my script here to use those blocks and pass them into the
    %%%%kinematics model.
    %%%%PSEUDOCODE:
    %%%%    1.) Calc. ambient air properties at current altitude.
    %%%%    2.) Calc. drag based on current conditions.
    %%%%    3.) Pass drag coefficient, etc. into dynamics/kinematics model.
    %%%% ------------------------------------------------------------------
    [rho(i), T] = air_prop_calc(s(2,i)+launch_alt);   %get air properties at given altitude
    P_a(i) = rho(i)*287*T;        %Update current atmospheric P (Pa).
    v_norm(i) = norm(v(:, i));       %get previous velocity magnitude
    c = sqrt(1.4*287*T);     %calculate speed of sound at given temp
    M(i) = v_norm(i)/c;                 %calculate Mach number
    C_D(i) = Coefficient_of_Drag_Defiance(Drag_Model, s(2,i)+launch_alt, M(i), is_burning);          %get coefficient of drag
    S = pi/4*diam^2;	%Calc. cross-sectional area of rocket.
    F_D(i) = rho(i)/2*(v(1, i)^2 + v(2, i)^2)*S*C_D(i);    %Drag mag..
    P_dyn(i) = (rho(i)/2)*(v_norm(i))^2;                   %Dynamic Pressure 
    %Check if rocket is still on launch rail and change angle accordingly:
    if s(2, i) > rail_len
            theta = atan2d(v(1, i), v(2, i));
    end
    
    a(:, i) = (F(i) - F_D(i))*[sind(theta); cosd(theta)] ...
        /(m + m_f(i) + m_ox(i) + m_pr) - [0; g0];          %Calc. accel.
    v(:, i+1) = v(:, i) + a(:, i)*dt;               %Calc. new velocity.
    s(:, i+1) = s(:, i) + v(:, i)*dt + a(:, i)/2*dt^2;	%Get new position.

    Total_Impulse = Total_Impulse + F(i)*dt;
    
    a_norm(i) = norm(a(:,i));
    compressive_force(i) = F(i) + F_D(i) + g0*(m + m_f(i) + m_ox(i) + m_pr)*cosd(theta); %Compressive load on the vehicle
    % Thrust, drag, gravity though center axis; 
    v_D = v_D + F_D(i)*cosd(theta)/(m + m_f(i) + m_ox(i)+m_pr)*dt;   %Accum. drag loss.
	%%%% ------------------------------------------------------------------
    
    i = i + 1;      %Increment counter index.
end


% Calculate m_dot_N2 using a central difference scheme
m_dot_pr = zeros(its, 1);
m_dot_pr(1, 1) = MW_N2*(n_pr_store(2, 1) - n_pr_store(1, 1))/dt;

for j = 2:burn_step-1
    m_dot_pr(j, 1) = MW_N2*(n_pr_store(j+1, 1) - n_pr_store(j-1, 1))/(2*dt);
end

m_dot_pr(burn_step, 1) = MW_N2*(n_pr_store(burn_step, 1) - n_pr_store(j, 1))/dt;


%Trim monitored variables:
m_dot_ox = m_dot_ox(1:i);       %Ox mass flow rate
F = F(1:i);                     %Thrust
Pcc = Pcc(1:i);                 %Combustion Chamber Pressure
P_ox = P_ox(1:i);               %Ox tank pressure
T_ox = T_ox(1:i);               %Ox tank temperature
Tcc = Tcc(1:i);                 %Combustion chamber temperature
OF = OF(1:i);                   %Ox/Fuel mass ratio
C_D = C_D(1:i-1);               %Drag Coef
F_D = F_D(1:i);                 %Drag force
a = a(:, 1:i);                  %Acceleration components
v = v(:, 1:i);                  %Velocity components
v_norm = v_norm(1:i-1);         %Velocity magnitude
s = s(:, 1:i);                  %Position components
m_f = m_f(1:i);                 %Fuel mass
m_ox = m_ox(1:i);               %ox mass
M = M(1:i-1);                   %Mach Number
t = (0:dt:(i - 1)*dt);          %Time vector
G_ox = G_ox(1:i);               %Oxidizer Mass Flux
rho = rho(1:i);                 %Air Density
P_a = P_a(1:i);                 %Atmospheric Pressure
compressive_force = compressive_force(1:i); %Compressive load on the rocket due to acceleration and aerodynamic pressure
a_norm = [sqrt(sum(a(:,1:burn_step).^2)),-sqrt(sum(a(:,burn_step+1:end).^2))];    %Magnitude of Acceleration
n_pr_store = n_pr_store(1:i);   %Pressurant amount
m_dot_pr = m_dot_pr(1:i);       %Pressurant mass flow rate


%% Plot Outputs/Monitored Variables %%
close all

figure(1)
subplot(2,1,1);
plot(s(2,1:end-1), rho(1:end-1),'linewidth',2);
xlabel('Altitude (m)')
ylabel('Air Density (kg/m^3)')
grid minor

subplot(2,1,2)
plot(t(1:end-1), P_a(1:end-1),'linewidth',2);
xlabel('time (s)')
ylabel('Atmospheric Pressure (Pa)')
grid minor

figure(2)
subplot(3,1,1);
plot(t(1:end-1),F_D(1:end-1),...
    'linewidth',2);
xlabel('time (s)')
ylabel('Force of Drag (N)')
grid minor

subplot(3,1,2);
plot(t(1:end-1),C_D,...
    'linewidth',2);
xlabel('time (s)')
ylabel('Drag Coefficient')
grid minor

subplot(3,1,3);
plot(M(1:end-1),C_D(1:end-1),...
    'linewidth',2);
xlabel('Mach Number (s)')
ylabel('Drag Coefficient')
grid minor

figure(3)
subplot(2,2,1);
plot(t(1:2*burn_step),m_dot_ox(1:2*burn_step),...
    'linewidth',2);
xlabel('time (s)')
ylabel('Rate of m_{ox} (kg/s)')
xlim([0,2*burn_step*dt])
grid minor

subplot(2,2,2);
plot(t(1:2*burn_step),G_ox(1:2*burn_step),...
    'linewidth',2);
xlabel('time (s)')
ylabel('Oxidizer Mass Flux (kg/s)')
xlim([0,2*burn_step*dt])
grid minor

subplot(2,2,3);
plot(t(1:2*burn_step),OF(1:2*burn_step),...
    'linewidth',2);
xlabel('time (s)')
ylabel('O/F ratio')
xlim([0,2*burn_step*dt])
grid minor

subplot(2,2,4);
plot(t(1:2*burn_step),Tcc(1:2*burn_step),...
    'linewidth',2);
xlabel('time (s)')
ylabel('Flame Temperature (K)')
xlim([0,2*burn_step*dt])
grid minor

figure(4)
subplot(2,1,1);
plot(t(1:2*burn_step),F(1:2*burn_step),...
    'linewidth',2);
ylabel('Thrust (N)')
xlabel('time (s)')
grid minor

subplot(2,1,2);
plot(t(1:2*burn_step),P_ox(1:2*burn_step),t(1:2*burn_step),Pcc(1:2*burn_step),...
    'linewidth',2);
legend('Ox Tank Pressure','Combustion Chamber Pressure')
xlabel('time (s)')
ylabel('Pressure (Pa)')
grid minor

figure(5)
subplot(2,1,1);
plot(t(1:end-1), v_norm,...
    'linewidth',2);
xlabel('time (s)')
ylabel('Axial Velocity (m/s)')
grid minor

subplot(2,1,2);
plot(t(1:end-1), a_norm(1:end-1),...
    'linewidth',2);
xlabel('time (s)')
ylabel('Axial Acceleration (m/s^2)')
grid minor

figure(6)
plot(t, s(1, :), t, s(2, :),...
    'linewidth',2);
legend('Downrange Distance','Altitude', 'Location','northwest')
xlabel('time (s)')
ylabel('Distance (m)')
grid minor

disp('Results')
disp(' ')
disp(['Burn Time: ' num2str(t(burn_step)) ' sec'])
disp(['Apogee: ' num2str(s(2,end)) ' m'])
disp(['Downrange Distance at Apogee: ' num2str(s(1,end)) ' m'])
disp(['Burnout Altitude: ' num2str(burnout_alt) ' m'])
disp(['Percent of fuel remaining: ' num2str(m_f(burn_step)/m_f(1)) ' %'])
disp(['Mass of nitrogen pressurant required: ' num2str(m_pr) ' kg'])
%disp(['Volume of nitrogen pressurant at 4500 psi: ' num2str(vol_pr) ' m^3'])
disp(' ')
disp(['TWR at launch: ' num2str(F(5)/(g0*(m + m_f(5) + m_ox(5))))])
disp(['Max Acceleration: ' num2str(max(a_norm(5:end))) ' m/sec^2'])
disp(['Max Velocity: ' num2str(max(v_norm)) ' m/sec'])
disp(['Max Mach: ' num2str(max(M))])
disp(['Horizontal Velocity at Apogee: ' num2str(v_norm(end)) ' m/sec'])
disp(['Total Impulse of the Vehicle: ' num2str(Total_Impulse) ' Ns'])
disp(' ')
disp(['Max Drag Force: ' num2str(max(F_D)) ' N'])
disp(['Total Drag Losses: ' num2str(v_D) ' m/s'])
disp(['Max Axial Load (i.e. acceleration/drag/gravity): ' num2str(max(compressive_force)) ' N '])
[maxQ, ind] = max(P_dyn);
disp(['Max Dynamic Pressure: ' num2str(maxQ) ' Pa at ' num2str(ind*dt) 's']);
