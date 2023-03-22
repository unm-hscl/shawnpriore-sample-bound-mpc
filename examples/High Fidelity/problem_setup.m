% add path
addpath './Helper Functions'
addpath './Data Collection'

%%%%%%%%%%%%%%%%
% get dynamics
%%%%%%%%%%%%%%%%

sampling_period = 60; % sec 
time_horizon      = 5;
[Ad_concat, Bd_concat, Cd_concat] = get_cwh(sampling_period, time_horizon);

x_size = size(Ad_concat,1);
u_size = size(Bd_concat,2);


%%%%%%%%%%%%%%%%
% relative positions
%%%%%%%%%%%%%%%%

sat_setup;
x_0 = eci2cwh(x_0_chief, x_0_dep_1);
x_mean_no_input = Ad_concat * x_0;


%%%%%%%%%%%%%%%%
% target set
% format: x, y, z, x., y., z.
%%%%%%%%%%%%%%%%

G_k = [-1,  0,  1, 0, 0, 0;
       -1,  1,  0, 0, 0, 0;
       -1,  0, -1, 0, 0, 0;
       -1, -1,  0, 0, 0, 0;
        1,  0,  0, 0, 0, 0];
h_k  = [0;0;0;0;10];

G_N = kron(eye(6), [1;-1]);
h_N = [2; 0; ones(4,1); 0.5 * ones(6,1)];

G = blkdiag( kron(eye(time_horizon-1), G_k), G_N);
h = [kron(ones(time_horizon-1,1),h_k); h_N];
                                          
% get number of half space constraints                      
n_lin_state = size(G,1);

%%%%%%%%%%%%%%%%
% admissable input sets
%%%%%%%%%%%%%%%%

u_max = 0.04;
input_space_A = kron(eye(time_horizon), [eye(3); -eye(3)]);
input_space_b = u_max * ones(time_horizon*6,1);

%%%%%%%%%%%%%%%%
% probabilistic constraints
%%%%%%%%%%%%%%%%

safety_target         = 0.05;  % in target set

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimizer choices
%%%%%%%%%%%%%%%%%%%%%%%%%%%

cvx_solver gurobi;
cvx_solver_settings('TimeLimit', 1800);
cvx_precision default;

