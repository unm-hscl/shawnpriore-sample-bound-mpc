% add path
addpath './Helper Functions'
addpath './Data Collection'

%%%%%%%%%%%%%%%%
% relative positions
%%%%%%%%%%%%%%%%
sat_setup;
x_0_a = eci2cwh(x_0_chief, x_0_dep_1);
x_0_b = eci2cwh(x_0_chief, x_0_dep_2);
x_0_c = eci2cwh(x_0_chief, x_0_dep_3);

%%%%%%%%%%%%%%%%
% get dynamics
%%%%%%%%%%%%%%%%
sampling_period = 60; % sec 
time_steps      = 5;
[A_concat, C_concat, D_concat] = get_cwh(sampling_period, time_steps);

x_size = size(A_concat,1);
u_size = size(C_concat,2);

%%%%%%%%%%%%%%%%
% terminal target sets
% format: x, y, z, x., y., z.
%%%%%%%%%%%%%%%%

vmax = 0.25;

target_set_a = Polyhedron('lb', [ -2.5;  2.5;  7.5; -vmax; -vmax; -vmax], ...
                          'ub', [  2.5;  7.5; 12.5;  vmax;  vmax;  vmax]);  
target_set_b = Polyhedron('lb', [ -7.5;  -10; -2.5; -vmax; -vmax; -vmax], ... 
                          'ub', [ -2.5;   -5;  2.5;  vmax;  vmax;  vmax]);  
target_set_c = Polyhedron('lb', [ -2.5; -2.5;-12.5; -vmax; -vmax; -vmax], ...
                          'ub', [  2.5;  2.5; -7.5;  vmax;  vmax;  vmax]);  
                      
n_lin_state = size(target_set_c.A,1);

%%%%%%%%%%%%%%%%
% admissable input sets
% format: U_x, U_y, U_z
%%%%%%%%%%%%%%%%
u_max = 4;
input_space = Polyhedron('lb', [-u_max; -u_max; -u_max], ... 
                         'ub', [ u_max;  u_max;  u_max]);                         

input_space_A = kron( eye(time_steps), input_space.A);
input_space_b = kron( ones(time_steps,1), input_space.b);

%%%%%%%%%%%%%%%%
% probabilistic constraints
%%%%%%%%%%%%%%%%
safety_target           = 0.95;  % in target set
safety_collision_1_v    = 0.95;  % deputy and chief
safety_collision_2_v    = 0.95;  % deputy and deputy

r = 10; % collision avoid region radius
S = [eye(3), zeros(3)];
Sk = zeros(3, size(A_concat,2), time_steps);
for i = 1:time_steps
    Sk(:, 6*(i-1) + (1:6), i) = S;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sample data and statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_a = readmatrix('./Data Collection/a_dist.dat');
data_b = readmatrix('./Data Collection/b_dist.dat');
data_c = readmatrix('./Data Collection/c_dist.dat');

samples = size(data_a,2);

[W_a_mean, W_a_var] = get_sample_stats_1v(data_a);
[W_b_mean, W_b_var] = get_sample_stats_1v(data_b);
[W_c_mean, W_c_var] = get_sample_stats_1v(data_c);

[means_ab, Chol_vars_ab] = get_sample_stats_2v(S, D_concat, data_a, data_b, time_steps);
[means_ac, Chol_vars_ac] = get_sample_stats_2v(S, D_concat, data_a, data_c, time_steps);
[means_bc, Chol_vars_bc] = get_sample_stats_2v(S, D_concat, data_b, data_c, time_steps);

[means_a, Chol_vars_a] = get_sample_stats_2v(S, D_concat, data_a, zeros(size(data_a,1), size(data_a,2)), time_steps);
[means_b, Chol_vars_b] = get_sample_stats_2v(S, D_concat, data_b, zeros(size(data_b,1), size(data_b,2)), time_steps);
[means_c, Chol_vars_c] = get_sample_stats_2v(S, D_concat, data_c, zeros(size(data_c,1), size(data_c,2)), time_steps);

W_a_max = D_concat* max(data_a, [], 2);
W_b_max = D_concat* max(data_b, [], 2);
W_c_max = D_concat* max(data_c, [], 2);
W_a_min = D_concat* min(data_a, [], 2);
W_b_min = D_concat* min(data_b, [], 2);
W_c_min = D_concat* min(data_c, [], 2);

W_ab_max = D_concat* max(data_a-data_b, [], 2);
W_ac_max = D_concat* max(data_a-data_c, [], 2);
W_bc_max = D_concat* max(data_b-data_c, [], 2);
W_ab_min = D_concat* min(data_a-data_b, [], 2);
W_ac_min = D_concat* min(data_a-data_c, [], 2);
W_bc_min = D_concat* min(data_b-data_c, [], 2);


clear data_*;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimizer choices
%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter_max = 100;
tau = 1.1.^(1:iter_max);
epsilon_lambda  = 1e-8;
epsilon_dc      = 1e-5;
cvx_solver gurobi;
cvx_precision default;

