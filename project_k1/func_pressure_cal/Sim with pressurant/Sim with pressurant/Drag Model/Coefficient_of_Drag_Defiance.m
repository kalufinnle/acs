function [C_D] = Coefficient_of_Drag_Defiance(Drag_Model, h, M, burning)
%COEFFICIENT_OF_DRAG Interpolates to find current drag coefficient.
%   Returns drag coefficient of a rocket based on its current altitude and 
%   Mach # by interpolating the data specified by the input drag model.
%   INPUTS:
%       "Drag_Model" = ___________
%       "h" = scalar for current altitude [m].
%       "M" = scalar current Mach #.
%       "burning" = scalar logical where TRUE => engine is burning.
%   OUTPUT:
%       "C_D" = scalar drag coefficient.
% Written by Andreas Marquis and the propulsion gang - January 2017
% Revised: Saturday November 30, 2019 by Jacobious Maximus McGee
%------------------

if h > 19800
    h = 19800;   %temporary fix
end

if burning      %Reduce base drag while engine is burning.
    Cd_sea_level = interp1(Drag_Model(:,1),Drag_Model(:,3),M,'linear','extrap');
    Cd_65k = interp1(Drag_Model(:,1),Drag_Model(:,5),M,'linear','extrap');
    
    C_D = interp1([0; 19800],[Cd_sea_level; Cd_65k], h); %19800 m = 65 000 ft.
else            %Use normal base drag after burnout.
    Cd_sea_level = interp1(Drag_Model(:,1),Drag_Model(:,2),M,'linear','extrap');
    Cd_65k = interp1(Drag_Model(:,1),Drag_Model(:,4),M,'linear','extrap');
    
    C_D = interp1([0;19800],[Cd_sea_level;Cd_65k],h); %19800 m = 65 000 ft.
end

end
