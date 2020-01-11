function [m_dot_ox, Pcc, P_ox, Tcc, OF, m_f, m_ox, t, G_ox, F, i, burn_step] ...
    = engine_sim(imported)

%General parameters:
its = imported(1);          %Est'd # iterations.
dt = imported(2);           %Time step size (s).
comb_eff = imported(20);    %Combustive efficiency.
nozz_eff = imported(21);    %Nozzle efficiency.
Pa = imported(3);          %Ambient pressure at launch altitude (Pa).
T_a = imported(4);          %Ambient temperature at launch altitude (K).
rho_fuel = imported(6);     %Mass density of fuel (kg/m^3).

%CEA interpolation data:
num_rows = 356;     %# rows to import in CEA data file (= last row # - 1).
CEA_file = 'CEA_data.csv';
CEA = CEA_data(CEA_file,num_rows);

%Initialize variables and constants:
i = 1;              %Loop index reset.
g0 = 9.81;          %Acceleration due to gravity at sea-level (m/s^2).
R_u = 8314;         %Universal gas constant (J/kmol*K).
MW_ox = 44.013;     %Molecular weight/mass of N2O (kg/kmol).
f_D = 0.03;
L_fs = 0.00889;
d_fs = 0.0254;
A_inj = imported(18);
Cd = imported(19);
F = zeros(its, 1);  %Thrust (N).
Isp = F;            %Specific impulse (s).
OF = F;             %OF ratio.
Pcc = F;            %Combustion chamber pressure (Pa).
Pcc(1) = Pa;
Tcc = F;            %Combustion chamber temperature (K).
m_dot_ox = F;       %Ox mass flow rate (kg/s).
m_dot_f = F;        %Fuel mass flow rate (kg/s).
P_ox = F;           %Ox tank pressure array (Pa).
m_f = F;            %Fuel mass array (kg).
m_f(1) = imported(5);       %Initial fuel mass (kg).
I_total = 0;        %Monitor total impulse. 
G_ox = F;           %Oxidizer Flux

%Ox tank parameters:
m_ox = F;
m_ox(1) = imported(7);      %Loaded liquid ox mass in tank (kg).
T_ox = T_a;                 %Initial ox tank temperature (K).
V_tank = imported(8);      %Internal volume of ox tank (m^3).
m_tank = imported(9);      %Dry mass of ox tank (kg).
P_ox(1) = imported(22);     %Initial ox tank pressure (Pa). 4826330.105
Q = [2.781; 0.27244; 309.57; 0.2882];
        %Coefficients for molar specific vol. of liq. N2O (m^3/kmol).
Vhat_l = Q(2)^(1 + (1 - T_ox/Q(3))^Q(4))/Q(1);
        %Initial molar specific vol. of liq. N2O (m^3/kmol).
n(1) = m_ox(1)/MW_ox;       %# moles of liquid ox; n(2) is gaseous moles.
n_tot = n(1) + P_ox(1)*(V_tank - Vhat_l*n(1))/(R_u*T_ox);
n(2) = P_ox(1)*(V_tank - n_tot*Vhat_l)/(R_u*T_ox - P_ox(1)*Vhat_l);
        %NOTE: Assuming P_sat = P_ox(1).
        
%Fuel core parameters:
a0 = imported(10);          %Regression rate coefficient (m/s^2).
n_reg = imported(11);       %Regression rate exponent.
num_ports = imported(12);   %# combustion ports in fuel grain.
L = imported(13);           %Fuel core length (m).
w = imported(14);           %Initial fuel web thickness (m).
r_fo = imported(15);        %Inner combustion chamber radius (m).
A_cc = pi*r_fo^2;           %Combustion chamber cross-section area (m^2).

%Nozzle parameters:
A_t = imported(16);         %[Initial] throat area (m^2).
A_e = imported(17);         %Exit area (m^2).
e = @(Ma_e, k) (2/(k+1)*(1 + (k-1)/2*Ma_e^2))^((k+1)/(2*k-2))/Ma_e;
        %Function handle for nozzle area expansion ratio.
% e = @(Ma_cc, k) A_cc/A_e*Ma_cc/(1 + (k-1)/2*Ma_cc^2)^((k+1)/(2*k-2));
%         %Function handle for chamber-exit area-Mach relation.

%% Simulation Time Loop 
%Loop until engine stops burning
is_burning = w > 0 && Pcc(i) <= P_ox(i) && m_ox(i) > 0;
while is_burning
    is_burning = w > 0 && Pcc(i) <= P_ox(i) && m_ox(i) > 0;

    %If burn not finished, continue engine simulation:
    if is_burning
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
        
        F(i) = nozz_eff*(m_dot_prop*v_e + (P_e - Pa)*A_e);     %Thrust.
       
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
        Pcc(i) = Pa;
        P_ox(i) = Pa;
        m_f(i + 1) = m_f(i);
        m_ox(i + 1) = m_ox(i)*(m_f(i) > 0);
    end
    
    i = i + 1;      %Increment counter index.
end

%Trim monitored variables:
m_dot_ox = m_dot_ox(1:i);       %Ox mass flow rate
Pcc = Pcc(1:i);                 %Combustion Chamber Pressure
P_ox = P_ox(1:i);               %Ox tank pressure
Tcc = Tcc(1:i);                 %Combustion chamber temperature
OF = OF(1:i);                   %Ox/Fuel mass ratio
m_f = m_f(1:i);                 %Fuel mass
m_ox = m_ox(1:i);               %ox mass
t = (0:dt:(i - 1)*dt);          %Time vector
G_ox = G_ox(1:i);               %Oxidizer Mass Flux
F = F(1:i);                     %Thrust
end