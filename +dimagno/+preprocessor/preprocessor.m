%{
The preprocessor prepares the simulation data structure by running first
the default simrc file, then the simrc file provided by the user (if
any), and then overwriting with the input structure provided by the user
(if any).  

INPUT
* simrcfile: path to user simrc file to overwrite defaults
* userdata: data structure used to (partially) overwrite the data
  structure of the default smrc and the user simrc 

OUTPUT
* data: complete data structure with all the parameters of the problem
%}
function data = preprocessor(simrcfile,userdata)

%% Process input data
% Set defaults
data = dimagno.simrc; 

% simrcfile overwrites defaults
if exist('simrcfile','var') && ~isempty(simrcfile)
    [simrcfilepath,simrcfilename] = fileparts(simrcfile); 
    if ~isempty(simrcfilepath)
        currentdir = cd(simrcfilepath);
    end
    usersimrc = str2func(simrcfilename);
    if ~isempty(simrcfilepath)
        cd(currentdir);
    end
    data = usersimrc(data);
end

% userdata structure overwrites defaults 
if exist('userdata','var') && ~isempty(userdata)
    f = fieldnames(userdata); 
    for i = 1:length(f)
        data.(f{i}) = userdata.(f{i});
    end
end   

%% General
data.dimagno.datafile = fullfile(data.dimagno.simdir,'data.mat'); % path to file where full data structure will be saved
 
%% Logger
data.logger.logfile = fullfile(data.dimagno.simdir,'log.txt'); % path to file where the log will be saved
