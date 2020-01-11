%Main script for simulating full rocket performance. Uses engine simulation
%and 2D kinematics simulation to determine altitude of rocket. Employs
%Euler's Method to step variables through time. Uses pre-combustion values 
%in kinematics model to provide conservative estimate. Ignores ignition 
%transient and throat erosion. Plots monitored variables against time.
%Original scripts: Thomas S.H. Leung
%Revisions: Mitchell L. Passarelli
%Version Date: Saturday, March 11, 2017

%%%%%GENRAL NOTES:
%%%%%   In Thomas' file "Project14.m" I use full line comments to track
%%%%%   lines I have either accounted for or do not need in this file.

%%%%%NEW MEMBER TO-DO:
%%%%%   1.) Implement vapour burn.
%%%%%   2.) Implement cross-wind wind effects.

%%%%%SHORT-TERM TO-DO:
%%%%%   1.) See if I can use analytical expression for port diameter to
%%%%%   find burn time (by solving d_i - d(t) = 0) to get length of
%%%%%   matrices of monitored variables. Also see if I can couple it with
%%%%%   ox_tank model and pressure calculator to find the minumum of the
%%%%%   three burnout conditions.
%%%%%   2.) Test code by compaing with Thomas'. Make sure to check that the
%%%%%   results for the ox tank model are the same (esp. because of the use
%%%%%   of the n_l and n_g equations).

%%%%%LONG-TERM TO-DO (i.e., after finished):
%%%%%   1.) Look for a way to vectorize the loop (start with engine calcs).
%%%%%       Looks like only the ox tank model and fuel regression need to
%%%%%       be looped.
%%%%%   2.) Implement a check for choked flow and execute different
%%%%%       functions for choked & unchoked flow.

clear; clc; close all;

%% Import Inputs From Excel & Initialize Variables %%
%importing values into imported
imported = import_params('Preliminary Engine Design.xlsx', ...
    'Import Data_Deliverance II', 1, 26);

%General parameters:
its = imported(1);          %Est'd # iterations.
dt = imported(2);           %Time step size (s).
theta = imported(3);        %Launch angle, from vertical (deg).
rail_len = 10;              %Launch rail length (m).
comb_eff = imported(23);    %Combustive efficiency.
nozz_eff = imported(24);    %Nozzle efficiency.
Pa = imported(4);          %Ambient pressure at launch altitude (Pa).
T_a = imported(5)-8;          %Ambient temperature at launch altitude (K).
m = imported(6);            %Rocket dry mass (kg).
rho_fuel = imported(8);     %Mass density of fuel (kg/m^3).
launch_alt = imported(25);  %used to calculate atmospheric density and ...

%CEA interpolation data:
num_rows = 356;     %# rows to import in CEA data file (= last row # - 1).
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
g0 = 9.81;          %Acceleration due to gravity at sea-level (m/s^2).
R_u = 8314;         %Universal gas constant (J/kmol*K).
MW_ox = 44.013;     %Molecular weight/mass of N2O (kg/kmol).
v_D = 0;            %Drag [velocity] loss (m/s).
f_D = 0.03;
L_fs = 0.00889;
d_fs = 0.0254;
A_inj = imported(21);
Cd = imported(22);
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
C_D = F;            %Drag coefficient array
F_D = F;            %Drag force array (N).
m_f = F;            %Fuel mass array (kg).
M = F;              %Mach nnumber.
m_f(1) = imported(7);       %Initial fuel mass (kg).
diam = imported(20);        %Diameter of rocket (m).
dtheta = 0;         %change in the angle of rocket with respect to vertical 
                    %this is zero since we are not accoutning for crosswind yet    
I_total = 0;        %Monitor total impulse. 
v_norm = F;         %Velocity magnitude.
a_norm = F;         %Acceleration magnitude.
G_ox = F;           %Oxidizer Flux
rho = F;            %Air Density
P_a = F;            %Atmospheric Pressure
P_a(1) = Pa;
P_dyn = F;          %Dynamic Pressure

%Ox tank parameters:
m_ox = F;
m_ox(1) = imported(9);      %Loaded liquid ox mass in tank (kg).
T_ox = T_a;                 %Initial ox tank temperature (K).
V_tank = imported(10);      %Internal volume of ox tank (m^3).
m_tank = imported(11);      %Dry mass of ox tank (kg).
P_ox(1) = imported(26);     %Initial ox tank pressure (Pa). 4826330.105
Q = [2.781; 0.27244; 309.57; 0.2882];
        %Coefficients for molar specific vol. of liq. N2O (m^3/kmol).
Vhat_l = Q(2)^(1 + (1 - T_ox/Q(3))^Q(4))/Q(1);
        %Initial molar specific vol. of liq. N2O (m^3/kmol).
n(1) = m_ox(1)/MW_ox;       %# moles of liquid ox; n(2) is gaseous moles.
n_tot = n(1) + P_ox(1)*(V_tank - Vhat_l*n(1))/(R_u*T_ox);
n(2) = P_ox(1)*(V_tank - n_tot*Vhat_l)/(R_u*T_ox - P_ox(1)*Vhat_l);
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

addpath('Drag Model')
Drag_Model = Deliverance_AeroLab_Drag_Model();

%% Determine whether flow sys or injector plate primarily constricts flow:
%%%%%NOTE: Don't worry about this one yet. Realistically, the injector
%%%%%plate should always constrict the flow.


%% Simulation Time Loop (ignition to apogee) %%
%Loop until rocket reaches apogee (ignore condition on first iteration):
while s(2, i) > 0 || i == 1
    is_burning = w > 0 && Pcc(i) <= P_ox(i) && m_ox(i) > 0;
    
    %If burn not finished, continue engine simulation:
    if is_burning
        burnout_alt = s(2, i);    %Update burnout altitude.
        burn_step = i;            %Step when engine finishes its burn
        %Simulate blow-down in ox tank:
        if i == 1   %Uses initial guesses for initial Pcc and m_dot_ox.
            [m_dot_ox(i), dn, dT_ox, P_ox(i + 1)] = ox_tank_sim(...
                n, T_ox, V_tank, m_tank, 1.8, f_D, L_fs, d_fs, ...
                A_inj, Cd, 2.75e6);
        else        %Uses prev. iteration data.
            [m_dot_ox(i), dn, dT_ox, P_ox(i + 1)] = ox_tank_sim(...
                n, T_ox, V_tank, m_tank, m_dot_ox(i - 1), f_D, ...
                L_fs, d_fs, A_inj, Cd, Pcc(i - 1));
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
        if OF(i) > 7
            OF(i) = 7;
        end
        Tcc(i) = interp1(CEA(:,1), CEA(:,2), OF(i));
        gamma = interp1(CEA(:,1), CEA(:,4), OF(i));
        R = R_u/interp1(CEA(:,1), CEA(:,3), OF(i));
        char_vel = comb_eff*sqrt(gamma*R*Tcc(i))/gamma * ...
            ((gamma + 1)/2)^((gamma + 1)/(2*gamma - 2));
        %Eval. engine performance including nozzle:
        if i > 1   %Calc. combustion chamber P (avg to damp out osc.).
            Pcc(i) = (m_dot_prop*char_vel/A_t + Pcc(i - 1))/2;
        else
            Pcc(i) = m_dot_prop*char_vel/A_t/2;
        end
        Ma_e = fzero(@(Ma_e) A_e/A_t - e(Ma_e, gamma), 2.75);   %Exit Mach.
        P_e = Pcc(i)/(1 + (gamma-1)/2*Ma_e^2)^(gamma/(gamma-1));
                %Flow is always choked for Pcc >= 200 psi.
        
%         %Eval. engine performance when flow not choked:
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
        m_ox(i + 1) = m_ox(i) + dn(1)*MW_ox*dt; %Consume ox mass.
        n = n + dn*dt;                      %Consume moles of ox.
        T_ox = T_ox + dT_ox*dt;             %Step ox tank temp.
        %%%%Include updates to ox tank parameters here too.%%%%
 
    %Change here for tracking
    else
        Pcc(i) = P_a(i-1);
        P_ox(i) = P_a(i-1);
        m_f(i + 1) = m_f(i);
        m_ox(i + 1) = m_ox(i)*(m_f(i) > 0);
        n(2) = 0;
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
    C_D(i) = Coefficient_of_Drag_Deliverance(Drag_Model, s(2,i)+launch_alt, M(i), is_burning);          %get coefficient of drag
    S = pi/4*diam^2;	%Calc. cross-sectional area of rocket.
    F_D(i) = rho(i)/2*(v(1, i)^2 + v(2, i)^2)*S*C_D(i);    %Drag mag..
    P_dyn(i) = (rho(i)/2)*(v_norm(i))^2;                   %Dynamic Pressure
    
    %Check if rocket is still on launch rail and change angle accordingly:
    if s(2, i) > rail_len
            theta = atan2d(v(1, i), v(2, i));
    end
    
    a(:, i) = (F(i) - F_D(i))*[sind(theta); cosd(theta)] ...
        /(m + m_f(i) + m_ox(i) + n(2)*MW_ox) - [0; g0];          %Calc. accel.
    v(:, i+1) = v(:, i) + a(:, i)*dt;               %Calc. new velocity.
    s(:, i+1) = s(:, i) + v(:, i)*dt + a(:, i)/2*dt^2;	%Get new position.

    Total_Impulse = Total_Impulse + F(i)*dt;
    
    a_norm(i) = norm(a(:,i));
    % Thrust, drag, gravity though center axis; 
    v_D = v_D + F_D(i)*cosd(theta)/(m + m_f(i) + m_ox(i))*dt;   %Accum. drag loss.
	%%%% ------------------------------------------------------------------
    
    i = i + 1;      %Increment counter index.
end

%Trim monitored variables:
m_dot_ox = m_dot_ox(1:i);       %Ox mass flow rate
F = F(1:i);                     %Thrust
Pcc = Pcc(1:i);                 %Combustion Chamber Pressure
P_ox = P_ox(1:i);               %Ox tank pressure
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
a_norm = [sqrt(sum(a(:,1:burn_step).^2)),-sqrt(sum(a(:,burn_step+1:end).^2))];    %Magnitude of Acceleration
P_s = 0.5*rho(1:end-1).*v_norm.^2+ P_a(1:end-1);               %Stagnation Pressure

%% Plot Outputs/Monitored Variables %%
close all

figure(1)
plot(s(2,1:end-1), rho(1:end-1),'linewidth',2);
xlabel('Altitude (m)')
ylabel('Air Density (kg/m^3)')
grid minor

figure(2)
subplot(3,1,1);
plot(t(1:end-1), P_s,'linewidth',2);
xlabel('time (s)')
ylabel('Stagnation Pressure (Pa)')
grid minor

subplot(3,1,2)
plot(t(1:end-1), P_a(1:end-1),'linewidth',2);
xlabel('time (s)')
ylabel('Atmospheric Pressure (Pa)')
grid minor

subplot(3,1,3)
plot(t(1:end-1), M,'linewidth',2);
xlabel('time (s)')
ylabel('Mach Number')
grid minor

figure(3)
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

figure(4)
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

figure(5)
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

figure(6)
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

figure(7)
plot(t, s(1, :), t, s(2, :),...
    'linewidth',2);
legend('Downrange Distance','Altitude')
xlabel('time (s)')
ylabel('Distance (m)')
grid minor

disp('Results')
disp(' ')
disp(['Ballistic Flight Time: ' num2str(t(i)) ' sec'])
disp(' ')
disp(['Burn Time: ' num2str(t(burn_step)) ' sec'])
disp(['Apogee: ' num2str(max(s(2,:))) ' m'])
disp(['Downrange Distance at Apogee: ' num2str(s(1,end)) ' m'])
disp(['Burnout Altitude: ' num2str(burnout_alt) ' m'])
disp(['Percent of fuel remaining: ' num2str(m_f(burn_step)/m_f(1)) ' %'])
disp(' ')
disp(['TWR at launch: ' num2str(F(5)/(g0*(m + m_f(5) + m_ox(5))))])
disp(['Max Acceleration: ' num2str(max(a_norm(5:end))) ' m/sec^2'])
disp(['Max Velocity: ' num2str(max(v_norm)) ' m/sec'])
disp(['Greatest Mach: ' num2str(max(M))])
disp(['Total Impulse of the Vehicle: ' num2str(Total_Impulse) ' m/s'])
disp(' ')
disp(['Max Drag Force: ' num2str(max(F_D)) ' N'])
disp(['Total Drag Losses: ' num2str(v_D) ' m/s'])
[maxQ, ind] = max(P_dyn);
disp(['Max Dynamic Pressure: ' num2str(maxQ) ' Pa at ' num2str(ind*dt) 's']);