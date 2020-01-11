%--------------------------------------------------------------------------
%                   UTAT ROCKETRY: HOT FIRE TEST
%                  HYBRID ROCKET ENGINE SIMULATION
%--------------------------------------------------------------------------

%Main script for simulating rocket engine simulation, based on main
%rocket simulation originally developed by Thomas S.H. Leung, 
% Mitchell L. Passarelli, and Andreas Marquis.
%Original scripts: Thomas S.H. Leung Mitchell L. Passarelli, 
%and Andreas Marquis.
%Revision: Abdullah Azzam
%Version Date: Saturday, August 12, 2017

clear; clc; close all;

%% Import Inputs From Excel & Initialize Variables %%
%importing values into imported
imported = import_params('Preliminary Engine Design.xlsx', ...
    'Import Engine_Sim', 1, 27);
dt = imported(2);
%% Run Engine Sim
[m_dot_ox, Pcc, P_ox, Tcc, OF, m_f, m_ox, t, G_ox, F, i, burn_step] = ...
    engine_sim(imported);

%% Plot Outputs/Monitored Variables %%
close all

figure(1)
subplot(2,2,1);
plot(t,m_dot_ox,...
    'linewidth',2);
xlabel('time (s)')
ylabel('Rate of m_{ox} (kg/s)')
xlim([0,1.2*i*dt])
grid minor

subplot(2,2,2);
plot(t,G_ox,...
    'linewidth',2);
xlabel('time (s)')
ylabel('Oxidizer Mass Flux (kg/s)')
xlim([0,1.2*i*dt])
grid minor

subplot(2,2,3);
plot(t,OF,...
    'linewidth',2);
xlabel('time (s)')
ylabel('O/F ratio')
xlim([0,1.2*i*dt])
grid minor

subplot(2,2,4);
plot(t,Tcc,...
    'linewidth',2);
xlabel('time (s)')
ylabel('Flame Temperature (K)')
xlim([0,1.2*i*dt])
grid minor

figure(2)
subplot(2,1,1);
plot(t,F,...
    'linewidth',2);
ylabel('Thrust (N)')
xlabel('time (s)')
grid minor

subplot(2,1,2);
plot(t,P_ox,t,Pcc,...
    'linewidth',2);
legend('Ox Tank Pressure','Combustion Chamber Pressure')
xlabel('time (s)')
ylabel('Pressure (Pa)')
grid minor

disp('    Results    ')
disp('---------------')
disp(['Burn Time: ' num2str(t(burn_step)) ' sec'])
disp(' ')
disp(['Percent of fuel remaining: ' num2str(m_f(burn_step)/m_f(1)) ' %'])
