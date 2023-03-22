start_time = tic;
cvx_begin quiet

    variable U_s(u_size, 1); 
    
    minimize (U_s'*U_s)
    subject to
        for n_i = 1:samples
            G * (x_mean_no_input + Bd_concat * U_s + Cd_concat * data(:, n_i)) <= h;
        end
        input_space_A * U_s <= input_space_b;
        
cvx_end

time_scenario = toc(start_time);
Ex_dep_s = x_mean_no_input + Bd_concat * U_s;
p = verify(Ex_dep_s, Cd_concat, G, h);

fprintf('Scenario Method \n');
fprintf('%s ', cvx_status);
fprintf('with %i samples \n', samples);
fprintf('Optimal Value: %e \n', cvx_optval);
fprintf('Time to Solve: %f \n', time_scenario);
fprintf('Empirical Constraint Satisfaction: %f \n\n', p)