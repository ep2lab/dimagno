%{
Example user simrc file for the tests

INPUT AND OUTPUT
* data: parameter structure (optional as input; can be initially empty). 
%} 
function data = example_simrc(data)

%% General
data.dimagno.simdir = fullfile(pwd,'aaaa'); % directory where simulation files will be saved
