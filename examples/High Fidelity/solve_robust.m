%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial input guess
%%%%%%%%%%%%%%%%%%%%%%%%%%
Ua_p = zeros(u_size, 1);
Ub_p = zeros(u_size, 1);
Uc_p = zeros(u_size, 1);

Ex_A_r = Ex_A_no;
Ex_B_r = Ex_B_no;
Ex_C_r = Ex_C_no;

start_time = tic;
for iter = 1:iter_max 

    % update collision avoid 2 vehicle
    [norm_a_min, norm_grad_a_min] = update_g_1v_robust(Ex_A_r, S, C_concat, W_a_min, time_steps);
    [norm_b_min, norm_grad_b_min] = update_g_1v_robust(Ex_B_r, S, C_concat, W_b_min, time_steps);
    [norm_c_min, norm_grad_c_min] = update_g_1v_robust(Ex_C_r, S, C_concat, W_c_min, time_steps);
    [norm_a_max, norm_grad_a_max] = update_g_1v_robust(Ex_A_r, S, C_concat, W_a_max, time_steps);
    [norm_b_max, norm_grad_b_max] = update_g_1v_robust(Ex_B_r, S, C_concat, W_b_max, time_steps);
    [norm_c_max, norm_grad_c_max] = update_g_1v_robust(Ex_C_r, S, C_concat, W_c_max, time_steps);

    % update collision avoid 2 vehicle
    [norm_ab_min, norm_grad_ab_min] = update_g_2v_robust(Ex_A_r, Ex_B_r, S, C_concat, W_ab_min, time_steps);
    [norm_ab_max, norm_grad_ab_max] = update_g_2v_robust(Ex_A_r, Ex_B_r, S, C_concat, W_ab_max, time_steps);
    [norm_ac_min, norm_grad_ac_min] = update_g_2v_robust(Ex_A_r, Ex_C_r, S, C_concat, W_ac_min, time_steps);
    [norm_ac_max, norm_grad_ac_max] = update_g_2v_robust(Ex_A_r, Ex_C_r, S, C_concat, W_ac_max, time_steps);
    [norm_bc_min, norm_grad_bc_min] = update_g_2v_robust(Ex_B_r, Ex_C_r, S, C_concat, W_bc_min, time_steps);
    [norm_bc_max, norm_grad_bc_max] = update_g_2v_robust(Ex_B_r, Ex_C_r, S, C_concat, W_bc_max, time_steps);

    
    cvx_begin quiet
        variable Ex_A_r(x_size);
        variable Ex_B_r(x_size);
        variable Ex_C_r(x_size);
        
        variable Ua(u_size);
        variable Ub(u_size);
        variable Uc(u_size);
        
        % 2 vehicle collision avoid
        variable slack_ab_min(time_steps) nonnegative;
        variable slack_ac_min(time_steps) nonnegative;
        variable slack_bc_min(time_steps) nonnegative;
        variable slack_ab_max(time_steps) nonnegative;
        variable slack_ac_max(time_steps) nonnegative;
        variable slack_bc_max(time_steps) nonnegative;
        
        % 2 vehicle collision avoid
        variable slack_a_min(time_steps) nonnegative;
        variable slack_b_min(time_steps) nonnegative;
        variable slack_c_min(time_steps) nonnegative;
        variable slack_a_max(time_steps) nonnegative;
        variable slack_b_max(time_steps) nonnegative;
        variable slack_c_max(time_steps) nonnegative;
        
        % cost params
        variable sum_slack;
        variable quad_input_cost;

        minimize (tau(iter)*sum_slack + quad_input_cost)
        subject to
            %----------------------------
            % cost variables
            %----------------------------
            quad_input_cost >=  Ua' * Ua + Ub' * Ub + Uc' * Uc;

            sum_slack == sum(slack_ab_min + slack_ac_min + slack_bc_min + slack_ab_max + slack_ac_max + slack_bc_max + slack_a_min + slack_b_min + slack_c_min + slack_a_max + slack_b_max + slack_c_max ); 
                      
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
            
            for k = 1:time_steps
                
                %----------------------------
                % colission avoidance constraint (1 vehicle)
                %----------------------------
                 - norm_a_min(k) - norm_grad_a_min(k,:) * (Ua - Ua_p) - slack_a_min(k) <= -r^2;
                 - norm_a_max(k) - norm_grad_a_max(k,:) * (Ua - Ua_p) - slack_a_max(k) <= -r^2;
                 - norm_b_min(k) - norm_grad_b_min(k,:) * (Ub - Ub_p) - slack_b_min(k) <= -r^2;
                 - norm_b_max(k) - norm_grad_b_max(k,:) * (Ub - Ub_p) - slack_b_max(k) <= -r^2;
                 - norm_c_min(k) - norm_grad_c_min(k,:) * (Uc - Uc_p) - slack_c_min(k) <= -r^2;
                 - norm_c_max(k) - norm_grad_c_max(k,:) * (Uc - Uc_p) - slack_c_max(k) <= -r^2;
                
                %----------------------------
                % colission avoidance constraint (2 vehicle)
                %----------------------------
                - norm_ab_min(k) - norm_grad_ab_min(k,:) * [Ua - Ua_p; Ub - Ub_p] - slack_ab_min(k) + r^2 <= 0;
                - norm_ac_min(k) - norm_grad_ac_min(k,:) * [Ua - Ua_p; Uc - Uc_p] - slack_ac_min(k) + r^2 <= 0;
                - norm_bc_min(k) - norm_grad_bc_min(k,:) * [Ub - Ub_p; Uc - Uc_p] - slack_bc_min(k) + r^2 <= 0;
                - norm_ab_max(k) - norm_grad_ab_max(k,:) * [Ua - Ua_p; Ub - Ub_p] - slack_ab_max(k) + r^2 <= 0;
                - norm_ac_max(k) - norm_grad_ac_max(k,:) * [Ua - Ua_p; Uc - Uc_p] - slack_ac_max(k) + r^2 <= 0;
                - norm_bc_max(k) - norm_grad_bc_max(k,:) * [Ub - Ub_p; Uc - Uc_p] - slack_bc_max(k) + r^2 <= 0;
                
            end
            
            %----------------------------
            % terminal state constraint
            %----------------------------             
            target_set_a.A * (Ex_A_r(end-5:end) + W_a_min(end-5:end)) - target_set_a.b <= 0;
            target_set_a.A * (Ex_A_r(end-5:end) + W_a_max(end-5:end)) - target_set_a.b <= 0;
            target_set_b.A * (Ex_B_r(end-5:end) + W_b_min(end-5:end)) - target_set_b.b <= 0;
            target_set_b.A * (Ex_B_r(end-5:end) + W_b_max(end-5:end)) - target_set_b.b <= 0;
            target_set_c.A * (Ex_C_r(end-5:end) + W_c_min(end-5:end)) - target_set_c.b <= 0;
            target_set_c.A * (Ex_C_r(end-5:end) + W_c_max(end-5:end)) - target_set_c.b <= 0;
            
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