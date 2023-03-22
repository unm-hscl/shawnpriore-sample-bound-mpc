%% disturbance samples
big_M = 5000;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solve the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = tic;
cvx_begin quiet
    variable U_pc(u_size, 1);
    variable M_bin(samples, 1) binary;

    minimize (U_pc'*U_pc)
    subject to

        % u \in \mathcal{U_pc} 
        input_space_A * U_pc <= input_space_b;

        % mean in shrunk target set
        for i=1:samples
            G * (x_mean_no_input + Bd_concat * U_pc + Cd_concat * data(:, i)) - h <= M_bin(i) * big_M * ones(size(h,1),1);
        end
        sum(M_bin)/samples <= safety_target;
cvx_end


time_pc = toc(start_time);
Ex_dep_pc = x_mean_no_input + Bd_concat * U_pc;
p = verify(Ex_dep_pc, Cd_concat, G, h);

fprintf('Particle Method \n');
fprintf('%s ', cvx_status);
fprintf('with %i samples \n', samples);
fprintf('Optimal Value: %e \n', cvx_optval);
fprintf('Time to Solve: %f \n', time_pc);
fprintf('Empirical Constraint Satisfaction: %f \n\n', p)
