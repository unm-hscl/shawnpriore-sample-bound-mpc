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

clear data_*;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% power approx
%%%%%%%%%%%%%%%%%%%%%%%%%%
if (samples == 1000) && isfile('pow_func_1000.mat')
    load('pow_func_1000.mat');
else
    lam_min = sqrt(samples+1)*(cos(1/3*acos(-(samples-1)/(samples+1)))-0.5);
    pow_func = @(x) (sqrt(samples+1)+x).^2./(x.^2*samples+(sqrt(samples+1)+x).^2);
    [pow_func_m, pow_func_c] = function_affine(0, 1e-2, 1e4, lam_min, pow_func, 1e-3, lam_min);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% power approx
%%%%%%%%%%%%%%%%%%%%%%%%%%
hat_a = (1-safety_collision_2_v) / 15;
hat_lam = ((hat_a - 1) * sqrt(samples + 1) - sqrt(samples * hat_a * (1 - hat_a) * (samples + 1))) / (1 - hat_a * (samples + 1));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% holders
%%%%%%%%%%%%%%%%%%%%%%%%%%
input_cost_our_method = [1e10; zeros(iter_max,1)];
lambda_sum_our_method = [1e10; zeros(iter_max,1)];
total_cost_our_method = [1e20; zeros(iter_max,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial input guess
%%%%%%%%%%%%%%%%%%%%%%%%%%
Ua_p = zeros(u_size, 1);
Ub_p = zeros(u_size, 1);
Uc_p = zeros(u_size, 1);

Ex_A_no = A_concat * x_0_a;
Ex_B_no = A_concat * x_0_b;
Ex_C_no = A_concat * x_0_c;

Ex_A = Ex_A_no;
Ex_B = Ex_B_no;
Ex_C = Ex_C_no;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% precomputes
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ex[G * D(N) * W]
Ex_GDW_A = target_set_a.A * D_concat(end-5:end,:)*W_a_mean;
Ex_GDW_B = target_set_b.A * D_concat(end-5:end,:)*W_b_mean;
Ex_GDW_C = target_set_c.A * D_concat(end-5:end,:)*W_c_mean;

% Var[G * D(N) * W]
Var_GDW_A = sqrt(diag(target_set_a.A * D_concat(end-5:end,:) * W_a_var * D_concat(end-5:end,:)' * target_set_a.A'));
Var_GDW_B = sqrt(diag(target_set_b.A * D_concat(end-5:end,:) * W_b_var * D_concat(end-5:end,:)' * target_set_b.A'));
Var_GDW_C = sqrt(diag(target_set_c.A * D_concat(end-5:end,:) * W_c_var * D_concat(end-5:end,:)' * target_set_c.A'));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = tic;
for iter = 1:iter_max 

    % update collision avoid 2 vehicle
    [norm_a, norm_grad_a] = update_g_1v(Ex_A, S, C_concat, means_a, time_steps);
    [norm_b, norm_grad_b] = update_g_1v(Ex_B, S, C_concat, means_b, time_steps);
    [norm_c, norm_grad_c] = update_g_1v(Ex_C, S, C_concat, means_c, time_steps);

    % update collision avoid 2 vehicle
    [norm_ab, norm_grad_ab] = update_g_2v(Ex_A, Ex_B, S, C_concat, means_ab, time_steps);
    [norm_ac, norm_grad_ac] = update_g_2v(Ex_A, Ex_C, S, C_concat, means_ac, time_steps);
    [norm_bc, norm_grad_bc] = update_g_2v(Ex_B, Ex_C, S, C_concat, means_bc, time_steps);
    
    cvx_begin quiet
        variable Ex_A(x_size);
        variable Ex_B(x_size);
        variable Ex_C(x_size);
        
        variable Ua(u_size);
        variable Ub(u_size);
        variable Uc(u_size);
        
        % 2 vehicle collision avoid
        variable slack_ab(time_steps) nonnegative;
        variable slack_ac(time_steps) nonnegative;
        variable slack_bc(time_steps) nonnegative;
        
        % 2 vehicle collision avoid
        variable slack_a(time_steps) nonnegative;
        variable slack_b(time_steps) nonnegative;
        variable slack_c(time_steps) nonnegative;

        % targe set
        variable lambda_a(n_lin_state) nonnegative;
        variable lambda_b(n_lin_state) nonnegative;
        variable lambda_c(n_lin_state) nonnegative;
        variable prob_a(n_lin_state) nonnegative;
        variable prob_b(n_lin_state) nonnegative;
        variable prob_c(n_lin_state) nonnegative;
        
        % cost params
        variable sum_slack;
        variable quad_input_cost;

        minimize (tau(iter)*sum_slack + quad_input_cost)
        subject to
            %----------------------------
            % cost variables
            %----------------------------
            quad_input_cost >=  Ua' * Ua + Ub' * Ub + Uc' * Uc;

            sum_slack == sum(slack_a + slack_b + slack_c + slack_ab + slack_ac + slack_bc); 
                      
            
            %----------------------------
            % linear equations defining the state
            %----------------------------
            Ex_A == Ex_A_no + C_concat * Ua;
            Ex_B == Ex_B_no + C_concat * Ub;
            Ex_C == Ex_C_no + C_concat * Uc;

            %----------------------------
            % u \in \mathcal{U} 
            %----------------------------
            input_space_A * Ua <= input_space_b;
            input_space_A * Ub <= input_space_b;
            input_space_A * Uc <= input_space_b;
            
            
            
            slack_ab >= 0;
            slack_ac >= 0;
            slack_bc >= 0;
            slack_a >= 0;
            slack_b >= 0;
            slack_c >= 0;
            
            for k = 1:time_steps
                
                %----------------------------
                % colission avoidance constraint (1 vehicle)
                %----------------------------
                hat_lam * norm( Chol_vars_a(:,:,k) * [Sk(:,:,k) * Ex_A; 1] ) ...
                    - norm_a(k) ...
                    - norm_grad_a(k,:) * (Ua - Ua_p)...
                    - slack_a(k) + r^2 <= 0;
                hat_lam * norm( Chol_vars_b(:,:,k) * [Sk(:,:,k) * Ex_B; 1] ) ...
                    - norm_b(k) ...
                    - norm_grad_b(k,:) * (Ub - Ub_p)...
                    - slack_b(k) + r^2 <= 0;
                hat_lam * norm( Chol_vars_c(:,:,k) * [Sk(:,:,k) * Ex_C; 1] ) ...
                    - norm_c(k) ...
                    - norm_grad_c(k,:) * (Uc - Uc_p)...
                    - slack_c(k) + r^2 <= 0;
                
                %----------------------------
                % colission avoidance constraint (2 vehicle)
                %----------------------------
                hat_lam * norm( Chol_vars_ab(:,:,k) * [Sk(:,:,k) *(Ex_A - Ex_B); 1] ) ...
                    - norm_ab(k) ...
                    - norm_grad_ab(k,:) * [Ua - Ua_p; Ub - Ub_p] ...
                    - slack_ab(k) + r^2 <= 0;
                hat_lam * norm( Chol_vars_ac(:,:,k) * [Sk(:,:,k) *(Ex_A - Ex_C); 1] ) ...
                    - norm_ac(k) ...
                    - norm_grad_ac(k,:) * [Ua - Ua_p; Uc - Uc_p] ...
                    - slack_ac(k) + r^2 <= 0;
                hat_lam * norm( Chol_vars_bc(:,:,k) * [Sk(:,:,k) *(Ex_B - Ex_C); 1] ) ...
                    - norm_bc(k) ...
                    - norm_grad_bc(k,:) * [Ub - Ub_p; Uc - Uc_p] ...
                    - slack_bc(k) + r^2 <= 0;
                
            end
            
            %----------------------------
            % terminal state constraint
            %----------------------------
            
            lambda_a >= lam_min;
            lambda_b >= lam_min;
            lambda_c >= lam_min;
             
            % mean in shrunk target set
            target_set_a.A * Ex_A(end-5:end) + Ex_GDW_A + lambda_a .* Var_GDW_A - target_set_a.b <= 0;
            target_set_b.A * Ex_B(end-5:end) + Ex_GDW_B + lambda_b .* Var_GDW_B - target_set_b.b <= 0;
            target_set_c.A * Ex_C(end-5:end) + Ex_GDW_C + lambda_c .* Var_GDW_C - target_set_c.b <= 0;
            
            prob_a >= 0;
            prob_b >= 0;
            prob_c >= 0;
            for i = 1:(n_lin_state)
                prob_a(i) >= pow_func_m .* lambda_a(i) + pow_func_c;
                prob_b(i) >= pow_func_m .* lambda_b(i) + pow_func_c;
                prob_c(i) >= pow_func_m .* lambda_c(i) + pow_func_c;
            end
            sum(prob_a + prob_b + prob_c) <= 1-safety_target;
            
    cvx_end
    
    % update Costs
    input_cost_our_method(iter+1) = quad_input_cost;
    lambda_sum_our_method(iter+1) = sum_slack;
    total_cost_our_method(iter+1) = cvx_optval;

    % calculate convergence criteria
    conv_check = abs(quad_input_cost - input_cost_our_method(iter) + tau(iter)*(sum_slack - lambda_sum_our_method(iter)));

    % print statistics
    if ~ quiet
        fprintf('iteration: %d ', iter);
        fprintf('\t %f', cvx_optval);
        fprintf('\t %e', conv_check); 
        fprintf('\t %e', sum_slack);
        fprintf('\t %f', toc(start_time));
        fprintf('\t %s \n', cvx_status);
    end

    % check for solved status
    if strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved')
        % check for convergence
        if (conv_check <= epsilon_dc) && (sum_slack <= epsilon_lambda)                 
           break
        end

        % if not converged update previous answer to current answer
        Ua_p = Ua;
        Ub_p = Ub;
        Uc_p = Uc;

    % if results are NaN break before error
    elseif strcmpi(cvx_status, 'Failed') || strcmpi(cvx_status, 'Infeasible')
        break
    end
end
total_time_our_method = toc(start_time);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% print some useful information
%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n%s in ', cvx_status);
fprintf('%i itterations \n', iter);
fprintf('Computation time (sec): %f \n', total_time_our_method);
fprintf('Total Cost: %f \n', cvx_optval);
fprintf('Slack Cost: %f \n', sum_slack);
fprintf('Input Cost: %f \n', quad_input_cost);

if strcmpi(cvx_status, 'Failed') || strcmpi(cvx_status, 'Infeasible')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%
verify;
