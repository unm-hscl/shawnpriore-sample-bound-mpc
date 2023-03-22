%lam_min = sqrt(samples + 1) * ( cos( acos( (1-samples)/(samples+1) ) / 3) - 0.5);
lam_min = 1/sqrt(3);
pow_func = @(x) (sqrt(samples+1)+x).^2./(x.^2*samples+(sqrt(samples+1)+x).^2);
[pow_func_m, pow_func_c] = function_affine(0, 2.5e-2, 5e4, 4.8, pow_func, 1e-5, 4.8);

[E_w, Cov_w] = get_sample_stats_1v(data);

std_est = zeros(n_lin_state,1);
for n_i = 1:n_lin_state
    std_est(n_i) = sqrt( G(n_i, :) * Cd_concat * Cov_w * Cd_concat' * G(n_i, :)');
end

start_time = tic;
cvx_begin quiet

    variable U(u_size, 1); 
    variable lambda(n_lin_state, 1);
    variable lambda_prob(n_lin_state, 1);
    
    minimize (U'*U)
    subject to
    
        input_space_A * U <= input_space_b;
        G * (x_mean_no_input + Bd_concat * U + Cd_concat * E_w) + lambda .* std_est <= h;
        lambda >= lam_min;
        lambda_prob >= 0;
        for n_i = 1:(n_lin_state)
                lambda_prob(n_i) >= pow_func_m .* lambda(n_i) + pow_func_c;
        end
        sum(lambda_prob) <= safety_target;
cvx_end

time_proposed = toc(start_time);
Ex_dep = x_mean_no_input + Bd_concat * U;
p = verify(Ex_dep, Cd_concat, G, h);


fprintf('Proposed Method \n');
fprintf('%s ', cvx_status);
fprintf('with %i samples \n', samples);
fprintf('Optimal Value: %e \n', cvx_optval);
fprintf('Time to Solve: %f \n', time_proposed);
fprintf('Empirical Constraint Satisfaction: %f \n\n', p)