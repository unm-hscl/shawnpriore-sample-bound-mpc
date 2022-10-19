%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all
cvx_clear;

quiet = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem system
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% problem set up
time_horizon = 8;

[Ad_concat, Bd_concat, Wd_concat] = get_cwh(time_horizon);

% initial states
% format: x, y, z,  x., y., z.
x_0_b = [98;   1; -3; 0; 0; 0] ; % satellite A
x_0_a = [95;   0;  5; 0; 0; 0] ; % satellite B
x_0_c = [100; -5; -0.1; 0; 0; 0] ; % satellite C
x_0 = [x_0_a, x_0_b, x_0_c];
vehicles = 3;
combinations = vehicles * (vehicles-1) / 2;

% target sets
% format: x, y, z, x., y., z.

target_set_b = Polyhedron('lb', [ -2.5; -2.5;  7.5; -0.01; -0.01; -0.01], ...
                          'ub', [  2.5;  2.5; 12.5;  0.01;  0.01;  0.01]);  
target_set_a = Polyhedron('lb', [ -2.5; -2.5;-12.5; -0.01; -0.01; -0.01], ...
                          'ub', [  2.5;  2.5; -7.5;  0.01;  0.01;  0.01]);  
target_set_c = Polyhedron('lb', [-7.5;     5; -2.5; -0.01; -0.01; -0.01], ... 
                          'ub', [-2.5;    10;  2.5;  0.01;  0.01;  0.01]);   
                      
target_set_all_A = target_set_a.A;
target_set_all_b = [target_set_a.b, target_set_b.b, target_set_c.b];
                      
% get number of half space constraints                      
n_lin_state = size(target_set_all_A,1);
                      
% Input space
u_max = 3;
input_space = Polyhedron('lb', [-u_max; -u_max; -u_max], ... 
                         'ub', [ u_max;  u_max;  u_max]);                         

input_space_A = blkdiag(input_space.A);
for i=1:(time_horizon-1)
    input_space_A = blkdiag(input_space_A, input_space.A);
end

input_space_b = repmat(input_space.b, time_horizon,1);

% collision avoid region radius
r = 14;

% safety threshold
safety_target           = 0.8;  % in target set
safety_collision_2_v    = 0.8;  % intersatellite

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sample data and statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 500;
W_a = get_data('./data/a_dist.dat', c);
W_b = get_data('./data/b_dist.dat', c);
W_c = get_data('./data/c_dist.dat', c);

[W_a_mean, W_a_var] = get_sample_stats_1v(W_a);
[W_b_mean, W_b_var] = get_sample_stats_1v(W_b);
[W_c_mean, W_c_var] = get_sample_stats_1v(W_c);

P = [eye(3), zeros(3)];
[means_ab, vars_ab] = get_sample_stats_2v(Wd_concat*(W_a-W_b), P);
[means_ac, vars_ac] = get_sample_stats_2v(Wd_concat*(W_a-W_c), P);
[means_bc, vars_bc] = get_sample_stats_2v(Wd_concat*(W_b-W_c), P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve probelm
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cvx_solver gurobi
cvx_precision default


iter_max = 20;

% convergence perameters
epsilon_dc = 1e-2; % convergence in cost
epsilon_lambda = 1e-8; % convergence of sum of slack variables to zero

% cost of slack variable
tau_max = 1e6;
tau_mult = 5;
tau = min(tau_max * ones(iter_max,1),  tau_mult.^(0:(iter_max-1))');


solve_our_method


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make some plots
%%%%%%%%%%%%%%%%%%%%%%%%%%


colors = [224,   0,   0; % red
           30, 144  255; % dark blue
            0, 170,  85; % green
          118,   0, 168; % purple
           46,  52,  59; % grey
          236, 176,  31;  % yellow
           76, 189, 237; % light blue
          161,  19,  46  % dark red
           ] ./ 255;
       
plot_symbols = ['o', 'd', '^', 'h', 'v', '>', 'p', 's'];

F_xy = [1,2,3,4];

fig = figure();
fig.WindowState = 'maximized';

subplot(2,1,1)

hold on 

for i = 1:3
    plot([x_0(1,i);x_mean_our_method(1:6:end,i)], [x_0(2,i);x_mean_our_method(2:6:end,i)], '-', 'Color', colors(i,:), 'Marker', plot_symbols(i));
end

patch('Faces',F_xy,'Vertices', Polyhedron(target_set_a.A([1;2;7;8],1:2), target_set_a.b([1;2;7;8])).V,...
    'FaceColor', colors(1,:), ...
    'FaceAlpha', 0.1); 
patch('Faces',F_xy,'Vertices', Polyhedron(target_set_b.A([1;2;7;8],1:2), target_set_b.b([1;2;7;8])).V,...
    'FaceColor', colors(2,:), ...
    'FaceAlpha', 0.1); 
patch('Faces',F_xy,'Vertices', Polyhedron(target_set_c.A([1;2;7;8],1:2), target_set_c.b([1;2;7;8])).V,...
    'FaceColor', colors(3,:), ...
    'FaceAlpha', 0.1); 
xlabel('x')
ylabel('y')
hold off

subplot(2,1,2)

hold on 

for i = 1:3
    plot([x_0(1,i);x_mean_our_method(1:6:end,i)], [x_0(3,i);x_mean_our_method(3:6:end,i)], '-', 'Color', colors(i,:), 'Marker', plot_symbols(i));
end

patch('Faces',F_xy,'Vertices', Polyhedron(target_set_a.A([1;3;7;9],[1;3]), target_set_a.b([1;3;7;9])).V,...
    'FaceColor', colors(1,:), ...
    'FaceAlpha', 0.1); 
patch('Faces',F_xy,'Vertices', Polyhedron(target_set_b.A([1;3;7;9],[1;3]), target_set_b.b([1;3;7;9])).V,...
    'FaceColor', colors(2,:), ...
    'FaceAlpha', 0.1); 
patch('Faces',F_xy,'Vertices', Polyhedron(target_set_c.A([1;3;7;9],[1;3]), target_set_c.b([1;3;7;9])).V,...
    'FaceColor', colors(3,:), ...
    'FaceAlpha', 0.1); 
xlabel('x')
ylabel('z')
hold off


