addpath '../'
addpath '../Helper Functions/'

sat_setup;

sample_period   = 60; %sec
time_steps      = 5;
tspan_full = 0:sample_period:sample_period*time_steps;
tspan_simp = 0:sample_period:sample_period;

options_ode = odeset('AbsTol',1e-8,'RelTol',1e-6);

%rng(1); % training
rng(2); % testing
baseline_JD = 2459927; % Dec 13, 2022

samples = 1000;
samples_dep_1 = zeros(6*time_steps, samples);
samples_dep_2 = zeros(6*time_steps, samples);
samples_dep_3 = zeros(6*time_steps, samples);


for sample = 1:samples
    % pick new date
    JD = baseline_JD + unifrnd(-30,30);
    
    [~, x_chief_full] = ode45(@(t,x)full_eom(t, x, JD), tspan_full, x_0_chief, options_ode);
    [~, x_dep_1_full] = ode45(@(t,x)full_eom(t, x, JD), tspan_full, x_0_dep_1, options_ode);
    [~, x_dep_2_full] = ode45(@(t,x)full_eom(t, x, JD), tspan_full, x_0_dep_2, options_ode);
    [~, x_dep_3_full] = ode45(@(t,x)full_eom(t, x, JD), tspan_full, x_0_dep_3, options_ode);


    for i = 1:time_steps
        [~, x_chief_simp] = ode45(@(t,x)simple_eom(t, x), tspan_simp, x_chief_full(i,:), options_ode);
        [~, x_dep_1_simp] = ode45(@(t,x)simple_eom(t, x), tspan_simp, x_dep_1_full(i,:), options_ode);
        [~, x_dep_2_simp] = ode45(@(t,x)simple_eom(t, x), tspan_simp, x_dep_2_full(i,:), options_ode);
        [~, x_dep_3_simp] = ode45(@(t,x)simple_eom(t, x), tspan_simp, x_dep_3_full(i,:), options_ode);


        samples_dep_1((i-1)*6 + [1:6], sample) = eci2cwh(x_chief_full(i+1,:)', x_dep_1_full(i+1,:)') - eci2cwh(x_chief_simp(end, :)', x_dep_1_simp(end, :)');
        samples_dep_2((i-1)*6 + [1:6], sample) = eci2cwh(x_chief_full(i+1,:)', x_dep_2_full(i+1,:)') - eci2cwh(x_chief_simp(end, :)', x_dep_2_simp(end, :)');
        samples_dep_3((i-1)*6 + [1:6], sample) = eci2cwh(x_chief_full(i+1,:)', x_dep_3_full(i+1,:)') - eci2cwh(x_chief_simp(end, :)', x_dep_3_simp(end, :)');
    end
end

% writematrix(samples_dep_1, 'a_dist.dat');
% writematrix(samples_dep_2, 'b_dist.dat'); % training
% writematrix(samples_dep_3, 'c_dist.dat');
writematrix(samples_dep_1, 'a_dist_test.dat');
writematrix(samples_dep_2, 'b_dist_test.dat'); % testing
writematrix(samples_dep_3, 'c_dist_test.dat');
