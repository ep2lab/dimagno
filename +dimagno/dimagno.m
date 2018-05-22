%{
This is the main dimagno function, intended to be run by the user. It
essentially chains together the preprocessor, the solver, and then any
postprocessing functions, and takes care of logging and saving to 
disk.  

INPUT
* simrcfile: path to user simrc file to overwrite defaults
* userdata: data structure used to (partially) overwrite the data
  structure of the default smrc and the user simrc 

OUTPUT
* data: complete data structure with all the parameters of the problem
* solution: structure with the solution of the simulation 
%}
function [data,solution] = dimagno(simrcfile,userdata)

%% Code version (must update manually when codebase changes)
VERSION = '20180521';

%% Run data preprocessor to prepare data structure
if ~exist('simrcfile','var')
    simrcfile = [];
end
if ~exist('userdata','var')
    userdata = [];
end
data  = dimagno.preprocessor.preprocessor(simrcfile,userdata);

%% Create dir and save data
[~,~,~] = mkdir(data.dimagno.simdir); % [~,~,~] avoids warnings if dir exists 
save(data.dimagno.datafile,'data');

%% Startup messages
logger.title(['DIMAGNO version ',VERSION],10,data.logger);
logger.log(['Simulation directory: ',data.dimagno.simdir],'INF',5,data.logger);
logger.log(['Data saved to ',data.dimagno.datafile,' successfully.'],'INF',5,data.logger);
 
%% Main simulation loop
i = data.initialfront.i;
front = data.initialfront.front;
[~,~,~] = mkdir([data.dimagno.simdir,'/fronts']); % [~,~,~] avoids warnings if dir exists 
save(fullfile(data.dimagno.simdir,'fronts',[num2str(i),'.mat']),'front');
if data.dimagno.plotfront    
    figure;
    front.plot_front;
    drawnow;
end
logger.log(['Front ',num2str(i),' copied from initial front and saved to disk'],'INF',5,data.logger);
while true
    if dimagno.exit.(data.exit.exit)(data,front,i)            
        logger.log(['Exit condition for main simulation loop satisfied for front: ',num2str(i)],'INF',5,data.logger);
        break;
    end
    try % Advance front        
        i = i+1;
        logger.log(['Beginning integration of front: ',num2str(i)],'INF',5,data.logger);
        front = dimagno.solver.(data.solver.solver)(data,front);                
        save(fullfile(data.dimagno.simdir,'fronts',[num2str(i),'.mat']),'front');
        if data.dimagno.plotfront
           front.plot_front;    
           drawnow;
        end
        logger.log(['Front ',num2str(i),' calculated successfully and saved to disk'],'INF',5,data.logger);
    catch % Advance front routine failed        
        logger.log(['Front advance routine failed on front ',num2str(i),'. Exiting main iteration loop...'],'ERR',5,data.logger);
        break;
    end
end
 
%% Run postprocessor functions and save solution
postfunctions = data.postprocessor.postfunctions;  
solution = struct();
for ifn = 1:length(postfunctions)
    logger.log(['Running postprocessor function: ',postfunctions{ifn}],'INF',5,data.logger);        
    solution = dimagno.postprocessor.(postfunctions{ifn})(data,solution);
end
solutionfile = fullfile(data.dimagno.simdir,'post.mat');
save(solutionfile,'-struct','solution');    
logger.log(['Solution saved to ',solutionfile],'INF',5,data.logger); 

%% Farewell
logger.title('DIMAGNO execution finished.',10,data.logger);










