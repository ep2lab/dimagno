%{
This postprocessing function reads all fronts and creates nice arrays
with the simulation data.

MMM20180522
%}
% function solution = fronts_to_arrays(data,solution)

frontfiles = dir(fullfile(data.dimagno.simdir,'fronts','*.mat'));

for i = data.initialfront.i:length(frontfiles)    
    load(fullfile(data.dimagno.simdir,'fronts',frontfiles(i).name),'front'); % loads current front
    f = str2num(frontfiles(i).name(1:end-4)); % get actual number by stripping '.mat'
    solution.z(:,f) = front.z_;
    solution.r(:,f) = front.r_;
    solution.u(:,f) = front.u_;
    solution.v(:,f) = front.v_;
    solution.w(:,f) = front.w_;
    solution.phi(:,f) = front.phi_;
    solution.psi(:,f) = front.psi_;
    solution.Bz(:,f) = front.Bz_;
    solution.Br(:,f) = front.Br_;
    for j = 1:data.plasma.n_electrons
        solution.(['ne',num2str(j)])(:,f) = front.ne_{j};
        solution.(['Te',num2str(j)])(:,f) = front.Te_{j};
        solution.(['we',num2str(j)])(:,f) = front.we_{j};
        solution.(['etae',num2str(j)])(:,f) = front.etae_{j};
    end
    solution.c(:,f) = front.c_;    
end

%% Additional variables

% Mach number
solution.M(:,f) = front.M_;

% total jtheta 
solution.jtheta = solution.z*0; % allocate
for j=1:data.plasma.n_electrons
    solution.jtheta = solution.jtheta + solution.(['ne',num2str(j)]).*(data.plasma.electrons{j}.q*solution.(['we',num2str(j)]) + data.plasma.ions.q*solution.w);
end

