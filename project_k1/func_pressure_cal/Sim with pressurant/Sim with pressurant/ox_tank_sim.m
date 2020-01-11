function [m_dot, dn, dT, P] = ox_tank_sim(n, T, V, m_tank, m_dot_prev, ...
    f_D, L, d, A_inj, Cd, Pe, n_pr, pr_flg)

%OX_TANK_SIM Models thermochemical processes of blow-down in oxidizer tank.
%   Propellant Tank Pressurization Modeling for a Hybrid Rocket - Margaret
%   Mary Fernandez
%
%   INPUTS:
%       "n" = 1x2 vector of ox moles in tank (kmol).
%       "T" = ox tank temperature (K).
%       "V" = volume of ox tank (m^3).
%       "m_tank" = mass of oxidizer tank (kg).
%       "m_dot_prev" = ox mass flow rate of previous time step (kg/s).
%       "f_D" = Darcy friction factor of flow sys.
%       "L" = flow sys length (m).
%       "d" = flow sys hydraulic diameter (m).
%       "A_inj" = total orifice area of injector plate (m^2).
%       "Cd" = injector plate discharge coefficient.
%       "Pe" = pressure at exit of injector plate (Pa).
%       "n_pr" = amount of pressurant added to tank (kmol).
%       "pr_flg" = 0 for helium 1 for nitrogen.
%
%   OUTPUTS:
%       "m_dot" = scalar ox mass flow rate (kg/s).
%       "dn" = 1x2 differential change in ox moles (kmol).
%       "dT" = differential change in ox tank temperature (K).
%       "P" = scalar internal pressure of ox tank (Pa).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Perry's Chemical Engineers' Handbook Property Equations:
G1 = 96.512;           %Vapour pressure of N2O [Pa] coefficients
G2 = -4045;            % valid for Temp range [182.3 K - 309.57 K]
G3 = -12.277;
G4 = 2.886e-5;
G5 = 2;
Tc = 309.57;           % critical temperature of N2O [K]
J1 = 2.3215e7;         % heat of vaporization of N2O [J/kmol] coefficients
J2 = 0.384;            % valid for Temp range [182.3 K - 309.57 K]
J3 = 0;
J4 = 0;

% heat capacity of N2 at constant pressure [J/(kmol*K)] coefficients
% valid for Temp range [100 K - 1500 K]
C1 = pr_flg*0.28883e5 + (1-pr_flg)*0.2079e5; 
C2 = 0;              
C3 = 0;
C4 = 0;
C5 = 0;

D1 = 0.2934e5; % heat capacity of N2O gas at constant pressure [J/(kmol*K)] coefficients
D2 = 0.3236e5; % valid for Temp range [100 K - 1500 K]
D3 = 1.1238e3;
D4 = 0.2177e5;
D5 = 479.4;
E1 = 6.7556e4; % heat capacity of N2O liquid at constant pressure [J/(kmol*K)] coefficients
E2 = 5.4373e1; % valid for Temp range [182.3 K - 200 K]
E3 = 0;
E4 = 0;
E5 = 0;
Q1 = 2.781;    % molar specific volume of liquid N2O [m^3/kmol] coefficients
Q2 = 0.27244;
Q3 = 309.57;
Q4 = 0.2882;
R = 8314;      %Universal gas constant (J/kmol*K).
MW = 44.013;   %Molar mass of N2O (kg/kmol).

% Calculate thermodynamic parameters as functions of temperature:
Vhat_l      = Q2^(1 + (1 - T/Q3)^Q4)/Q1;                    %Molar specific vol. of liq. N2O [m^3/kmol]
CVhat_He    = C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4 - R;     %Molar c_V of He [J/(kmol*K)]
CVhat_g     = D1 + D2*((D3/T)/sinh(D3/T))^2 + ...
              D4*((D5/T)/cosh(D5/T))^2 - R;                 %Molar c_V of N2O gas [J/(kmol*K)]
CVhat_l     = E1 + E2*T + E3*T^2 + E4*T^3 + E5*T^4;         %Molar c_V of N2O liq.; approx. = c_P [J/(kmol*K)]
Tr          = T/Tc;                                         %Reduced temperature.
delta_Hv    = J1*(1 - Tr) ^ (J2 + J3*Tr + J4*Tr^2);         %Heat of vaporization of N2O [J/kmol]
P_sat       = exp(G1 + G2/T + G3*log(T) + G4*T^G5);         %Vapour P of N20 (Pa).
dP_sat      = (-G2/(T^2) + G3/T + G4*G5*T^(G5-1)) ...       %Derivative of vapour P with respect to T.
              * exp(G1 + G2/T + G3*log(T) + G4*T^G5);        
c_P         = (4.8 + 0.00322*T)*155.239;                    %Specific heat of tank, Aluminum [J/kg*K]

% Calc. pipe P losses in flow sys:
liq_density = MW/Vhat_l;                    % Density of liquid nitrous                               
A = pi*(d^2)/4;                             % Pipe area   
v_fs = m_dot_prev/(A*liq_density);          % Flow speed
deltaP = (f_D*L*liq_density*v_fs^2)/(2*d);  % Pressure loss

%Simplified expression definitions for solution:
P = (n_pr + n(2))*R*T/(V - n(1)*Vhat_l);
a = m_tank*c_P + n_pr*CVhat_He + n(2)*CVhat_g + n(1)*CVhat_l;
RT = P*Vhat_l;
e = -delta_Hv + R*T;
f = -Cd*A_inj*sqrt(abs(2/MW*(P - deltaP - Pe)/Vhat_l));
j = -Vhat_l*P_sat;
k = (V - n(1)*Vhat_l)*dP_sat;
l = R*T;
q = R*n(2);
Z = (-f*(-j*a + (q - k)*RT))/(a*(l + j) + (q - k)*(e - RT));
W = (-Z*(l*a + (q - k)*e))/(-j*a + (q - k)*RT);

%Time derivative of tank conditions
m_dot = -f*MW;          % Mass flow rate
dn    = [W, Z];         % Time derivative of nitrous amount
dT    = (RT*W + e*Z)/a;  % Time derivative of tank temperature

end
