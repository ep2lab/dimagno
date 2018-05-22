%{
Compendium of all tests for the function/class in the name of this file.

NOTES:
* The functions to be tested must be in the Matlab path. You can run the
  tests by executing 'runtests'. 
* Your working directory must be the
  directory where this test file is contained.  
%}
function tests = test_dimagno
    clc;
    tests = functiontests(localfunctions);     
end

%----------------------------------------------------------------------
%----------------------------------------------------------------------
%----------------------------------------------------------------------
 
function test_ic(~)
    field = magnetic_field.loop_2d;
    plasma = fluid_plasma.plasma;
    ic = dimagno.ic('plasma',plasma,'field',field);
    F = dimagno.front;
    F.r_ = linspace(0,1);
    F = ic.initial_front(F);
end 

%----------------------------------------------------------------------
function test_preprocessor_and_default_simrc(~)
    % Call simrc without arguments
    data = dimagno.simrc; 
    % Call simrc with arguments
    data.test = 'just a test string';
    data = dimagno.simrc(data);
    % Call preprocessor without arguments
    data = dimagno.preprocessor.preprocessor;
    % Call preprocessor with a user simrc file
    data = dimagno.preprocessor.preprocessor(fullfile('fixtures/example_simrc.m'));
    % Call preprocessor with a data structute
    data = dimagno.preprocessor.preprocessor('fixtures/example_simrc.m',struct('testfield',12345));
end

%----------------------------------------------------------------------

function test_dimagnomain(~)
    [data,solution] = dimagno.dimagno;
end

%----------------------------------------------------------------------

function test_postproc(~) % IMPORTANT! Default simulation must exist before running this test
    load('sims/data.mat');
    solution = dimagno.postprocessor.fronts_to_arrays(data,[]);    
    solution = dimagno.postprocessor.induced_field(data,solution);    
end
























