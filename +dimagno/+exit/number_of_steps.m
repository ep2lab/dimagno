function do_exit = number_of_steps(data,front,i)  
% returns 1 if i > sim_params.exit_condition_parameters.n_steps
do_exit = i >= data.exit.parameters.n_steps;
