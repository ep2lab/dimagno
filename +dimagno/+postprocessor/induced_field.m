%{
This postprocessing function computes the plasma-induced magnetic field
from the azimuthal currents in the plasma. It must be run after
fronts_to_arrays. 

MMM20180522
%}
function solution = induced_field(data,solution)
   
%% Create and prepare induced field object
solution.induced = magnetic_field.library_2d; % Field object (library_2d)
 
%% Interpolation grid
Xlim = data.postprocessor.postparameters.induced_field.Xlim;
Ylim = data.postprocessor.postparameters.induced_field.Ylim;
x_points = data.postprocessor.postparameters.induced_field.xpoints;
y_points = data.postprocessor.postparameters.induced_field.ypoints;

[solution.induced.ZZ,solution.induced.RR] = meshgrid(linspace(Xlim(1),Xlim(2),x_points),linspace(Ylim(1),Ylim(2),y_points));
 
% Interpolate total azimuthal current on the grid
I = scatteredInterpolant(solution.z(:),solution.r(:),solution.jtheta(:));
JTHETA = I(solution.induced.ZZ,solution.induced.RR);
JTHETA(data.applied.field_2d(solution.induced.ZZ,solution.induced.RR) > data.applied.field_2d(data.initialfront.front.z_(end),data.initialfront.front.r_(end)) ) = 0; % set to 0 the jtheta of points outside of the plasma.
clear I; % delete interpolant 
 
% Generate temporary loop_2d object
t = magnetic_field.loop_2d;
t.ZL = (solution.induced.ZZ(1,2) + solution.induced.ZZ(1,1))/2; % loop positioned at half a cell
t.I = 1/t.const.mu0; % THIS REMOVES mu0, so we only need to multiply the field times beta0 to be able to add it with the applied one 20130512
cell_area = (solution.induced.ZZ(1,2)-solution.induced.ZZ(1,1))*(solution.induced.RR(2,1)-solution.induced.RR(1,1));

solution.induced.PSI = solution.induced.ZZ.*0; % allocate
solution.induced.BZ = solution.induced.PSI;
solution.induced.BR = solution.induced.PSI;
for k = 1:y_points-1 % Run along grid column
    t.RL = (solution.induced.RR(k+1,1) + solution.induced.RR(k,1))/2; % Half a cell
    [psi,Bz,Br] = t.field_2d(solution.induced.ZZ(:,2:end),solution.induced.RR(:,2:end));
    for j = 1:x_points-1 % Run along current row
        IL = cell_area * (JTHETA(k,j)+JTHETA(k+1,j)+JTHETA(k,j+1)+JTHETA(k+1,j+1))/4; % interpolated current at loop
        if IL == 0 || isnan(IL)
            continue;   % skip 0 and NaN loops
        end
        solution.induced.PSI = solution.induced.PSI + IL*[psi(:,j:-1:1),psi(:,1:end+1-j)]; % slice before and after the current loop position (mirror) and add it
        solution.induced.BZ  = solution.induced.BZ  + IL*[ Bz(:,j:-1:1), Bz(:,1:end+1-j)];
        solution.induced.BR  = solution.induced.BR  + IL*[-Br(:,j:-1:1), Br(:,1:end+1-j)];
    end
end  
