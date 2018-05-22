%{
Creates default data structure containing the parameters of the problem.
Called by main program to set the defaults, before running user simrc
file. You may take this file as a reference list of defaults and as a
template that you can copy to create your own user simrc files. 

INPUT AND OUTPUT
* data: parameter structure (optional as input; can be initially empty). 
%} 
function data = simrc(data)

%% General
data.dimagno.simdir = fullfile(pwd,'sims'); % directory where simulation files will be saved

data.dimagno.plotfront = 1; % whether to plot current front position during simulation

data.dimagno.max_iter = 100; % maximum number of iterations
data.dimagno.tol = 1e-6; % Tolerance to exit iterations
data.dimagno.corrector = 1; % corrector in predictor-corrector on?

data.dimagno.B_on_ions = 1; % 1/0 whether to activate magnetic force on ions (if 0, B only affects electrons)
data.dimagno.beta0 = 0.1; % plasma beta at the origin, beta0 =  mu0 n T / Ba0^2, for iterations of the self-induced field
data.dimagno.beta_0_tolerance = 1e-4;

%% Logger 
data.logger.filedebuglevel = 3; % file debug level. A higher number prints less messages
data.logger.screendebuglevel = 3; % screen debug level. A higher number prints less messages
data.logger.linelength = 80; % maximum line length in the logs

%% Applied and total magnetic field
data.applied = magnetic_field.loop_2d; % Only the applied field
data.field = data.applied; % Total (applied plus plasma-induced) field

%% Plasma
data.plasma = fluid_plasma.plasma;

%% Initial conditions
data.ic = dimagno.ic('plasma',data.plasma,'field',data.field);

%% Initial front
data.initialfront.i = 1; % Index for the initial front
data.initialfront.front = dimagno.front('r_',linspace(0,data.ic.r_max,20).','z_',zeros(20,1));
    data.initialfront.front = data.ic.initial_front(data.initialfront.front);
      
%% Solver
data.solver.solver = 'dMoC_advance_front'; % Type of front advancer to use               
data.exit.exit = 'number_of_steps'; % Type of exit condition to use
data.exit.parameters = struct(... % Parameters for exit condition. 
            'n_steps', 20); 

%% Postprocessor
data.postprocessor.postfunctions = {'fronts_to_arrays','induced_field'}; % Cell array with the names of postprocessor functions to run after iteration process
data.postprocessor.postparameters.induced_field = struct(...
            'Xlim',[0,0.35],...
            'Ylim',[0,1.1],...
            'xpoints',20,...
            'ypoints',20);
        
