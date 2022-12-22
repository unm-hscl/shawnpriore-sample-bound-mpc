%%%%%%%%%%%%%%%%%%%%%%%%%%%
% expectations
%%%%%%%%%%%%%%%%%%%%%%%%%%%
normal_mean = zeros(6*time_steps,1);
normal_var = kron(eye(5), blkdiag(1e-5*eye(3),1e-8*eye(3)));
normal_var_chol = chol(normal_var);

means_1v = zeros(size(S,1)+1, size(S,1)+1, time_steps);
vars_1v = zeros(size(S,1)+1, size(S,1)+1, time_steps);

means_2v = zeros(size(S,1)+1, size(S,1)+1, time_steps);
vars_2v = zeros(size(S,1)+1, size(S,1)+1, time_steps);

for i = 1:time_steps
    index = (i-1)*6 + [1:6];
    Dk = D_concat(index, :);
    
    VAR_Zk = S * Dk * normal_var * Dk' * S';
    VAR_ZkTZk = 2 * trace(VAR_Zk^2);

    vars_1v(:,:,i) = chol(blkdiag(4*VAR_Zk, VAR_ZkTZk));
    vars_2v(:,:,i) = chol(blkdiag(8*VAR_Zk, 4*VAR_ZkTZk));

    EX_ZkTZk = trace(VAR_Zk);

    means_1v(:,:,i) = blkdiag(eye(size(S,1)), EX_ZkTZk);
    means_2v(:,:,i) = blkdiag(eye(size(S,1)), 2 * EX_ZkTZk);

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% power approx
%%%%%%%%%%%%%%%%%%%%%%%%%%

lam_min = 1/sqrt(3);
pow_func = @(x) 1./(x.^2 + 1);
[pow_func_m, pow_func_c] = function_affine(0, 1e-2, 1e4, lam_min, pow_func, 1e-3, lam_min);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% power approx
%%%%%%%%%%%%%%%%%%%%%%%%%%
hat_a = (1-safety_collision_2_v) / 15;
hat_lam = sqrt(1/hat_a -1);


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

Ex_A_r = Ex_A_no;
Ex_B_r = Ex_B_no;
Ex_C_r = Ex_C_no;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% precomputes
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ex[G * D(N) * W]
Ex_GDW_A = target_set_a.A * D_concat(end-5:end,:)*normal_mean;
Ex_GDW_B = target_set_b.A * D_concat(end-5:end,:)*normal_mean;
Ex_GDW_C = target_set_c.A * D_concat(end-5:end,:)*normal_mean;

% Var[G * D(N) * W]
Var_GDW_A = sqrt(diag(target_set_a.A * D_concat(end-5:end,:) * normal_var * D_concat(end-5:end,:)' * target_set_a.A'));
Var_GDW_B = sqrt(diag(target_set_b.A * D_concat(end-5:end,:) * normal_var * D_concat(end-5:end,:)' * target_set_b.A'));
Var_GDW_C = sqrt(diag(target_set_c.A * D_concat(end-5:end,:) * normal_var * D_concat(end-5:end,:)' * target_set_c.A'));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = tic;
for iter = 1:iter_max 

    % update collision avoid 2 vehicle
    [norm_a, norm_grad_a] = update_g_1v(Ex_A_r, S, C_concat, means_1v, time_steps);
    [norm_b, norm_grad_b] = update_g_1v(Ex_B_r, S, C_concat, means_1v, time_steps);
    [norm_c, norm_grad_c] = update_g_1v(Ex_C_r, S, C_concat, means_1v, time_steps);

    % update collision avoid 2 vehicle
    [norm_ab, norm_grad_ab] = update_g_2v(Ex_A_r, Ex_B_r, S, C_concat, means_2v, time_steps);
    [norm_ac, norm_grad_ac] = update_g_2v(Ex_A_r, Ex_C_r, S, C_concat, means_2v, time_steps);
    [norm_bc, norm_grad_bc] = update_g_2v(Ex_B_r, Ex_C_r, S, C_concat, means_2v, time_steps);
    
    cvx_begin quiet
        variable Ex_A_r(x_size);
        variable Ex_B_r(x_size);
        variable Ex_C_r(x_size);
        
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
            Ex_A_r == Ex_A_no + C_concat * Ua;
            Ex_B_r == Ex_B_no + C_concat * Ub;
            Ex_C_r == Ex_C_no + C_concat * Uc;

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
                hat_lam * norm( vars_1v(:,:,k) * [Sk(:,:,k) * Ex_A_r; 1] ) ...
                    - norm_a(k) ...
                    - norm_grad_a(k,:) * (Ua - Ua_p)...
                    - slack_a(k) + r^2 <= 0;
                hat_lam * norm( vars_1v(:,:,k) * [Sk(:,:,k) * Ex_B_r; 1] ) ...
                    - norm_b(k) ...
                    - norm_grad_b(k,:) * (Ub - Ub_p)...
                    - slack_b(k) + r^2 <= 0;
                hat_lam * norm( vars_1v(:,:,k) * [Sk(:,:,k) * Ex_C_r; 1] ) ...
                    - norm_c(k) ...
                    - norm_grad_c(k,:) * (Uc - Uc_p)...
                    - slack_c(k) + r^2 <= 0;
                
                %----------------------------
                % colission avoidance constraint (2 vehicle)
                %----------------------------
                hat_lam * norm( vars_2v(:,:,k) * [Sk(:,:,k) *(Ex_A_r - Ex_B_r); 1] ) ...
                    - norm_ab(k) ...
                    - norm_grad_ab(k,:) * [Ua - Ua_p; Ub - Ub_p] ...
                    - slack_ab(k) + r^2 <= 0;
                hat_lam * norm( vars_2v(:,:,k) * [Sk(:,:,k) *(Ex_A_r - Ex_C_r); 1] ) ...
                    - norm_ac(k) ...
                    - norm_grad_ac(k,:) * [Ua - Ua_p; Uc - Uc_p] ...
                    - slack_ac(k) + r^2 <= 0;
                hat_lam * norm( vars_2v(:,:,k) * [Sk(:,:,k) *(Ex_B_r - Ex_C_r); 1] ) ...
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
            target_set_a.A * Ex_A_r(end-5:end) + Ex_GDW_A + lambda_a .* Var_GDW_A - target_set_a.b <= 0;
            target_set_b.A * Ex_B_r(end-5:end) + Ex_GDW_B + lambda_b .* Var_GDW_B - target_set_b.b <= 0;
            target_set_c.A * Ex_C_r(end-5:end) + Ex_GDW_C + lambda_c .* Var_GDW_C - target_set_c.b <= 0;
            
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
