time_steps  = 5;
samples     = 1e4;

norm_var = chol(kron(eye(5), blkdiag(1e-5*eye(3),1e-8*eye(3))));

%rng(10); % training
rng(2); % testing

samples_dep_1 = (randn(samples, 6*time_steps) * norm_var)';
samples_dep_2 = (randn(samples, 6*time_steps) * norm_var)';
samples_dep_3 = (randn(samples, 6*time_steps) * norm_var)';


% writematrix(samples_dep_1, 'a_dist.dat');
% writematrix(samples_dep_2, 'b_dist.dat'); % training
% writematrix(samples_dep_3, 'c_dist.dat');
writematrix(samples_dep_1, 'a_dist_test.dat');
writematrix(samples_dep_2, 'b_dist_test.dat'); % testing
writematrix(samples_dep_3, 'c_dist_test.dat');
