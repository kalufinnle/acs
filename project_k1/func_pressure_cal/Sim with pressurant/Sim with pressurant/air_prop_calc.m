%This function calculates air density and temperature at a specified altitude.
%Author: Abdullah Azzam
%Date: 10/31/2016

function [air_den,air_temp] = air_prop_calc(h)
% Input is altitude in metres
% Output is air density kg/m3, temp in K

%Parameters to be used

P0=101325; %P0: standard level atmospheric pressure, 101325 Pa
T0 = 288.15;%T0: sea level standard temperature, 288.15 K
g=9.81; %g: gravitational acceleration, 9.81 m/s2
a1=-0.0065;%L: temperature lapse rate, 0.0065 K/m
Ru=8.31447;%R: universal gas constant, 8.31447 J/(mol.K)
M=0.0289644;%M: molar mass of dry air, 0.0289644 kg/mol
R = Ru/M;%R: Specific gas constant of air J/(kg.K)

%Temperature approximation inside troposphere at altitude h

if h < 11000 
    T = T0 + (a1*h);
    
    %Air density in kg/m3
    air_den = 1.225*(T/T0)^(-g/(a1*R)-1);
    
else %if h < 25000 
    T = T0 + (a1*11000);
    air_den0 = 1.225*((T0 + (a1*11000))/T0)^(-g/(a1*R)-1);
    
    air_den = air_den0*exp(-g/(R*T)*(h-11000));
end

air_temp = T;

end
