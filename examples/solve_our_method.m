%%%%%%%%%%%%%%%%%%%%%%%%%%
% power approx
%%%%%%%%%%%%%%%%%%%%%%%%%%
lam_min = sqrt(c+1)*(cos(1/3*acos(-(c-1)/(c+1)))-0.5);
pow_func = @(x) (sqrt(c+1)+x).^2./(x.^2*c+(sqrt(c+1)+x).^2);
[pow_func_m, pow_func_c] = function_affine(0, 1e-3, 1000, lam_min, pow_func, 1e-4, lam_min);

hat_a = (1-safety_collision_2_v) / 24;
hat_lam = ((hat_a-1)*sqrt(c+1)-sqrt(c*hat_a*(1-hat_a)*(c + 1)))/(1-hat_a*(c+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%
% holders
%%%%%%%%%%%%%%%%%%%%%%%%%%
% storage initial cost for convergence check
input_cost_our_method = [1e10; zeros(iter_max,1)];
lambda_sum_our_method = [1e10; zeros(iter_max,1)];
total_cost_our_method = [1e20; zeros(iter_max,1)];

% initial input guess
U_p = zeros(size(Bd_concat,2), 3);
x_mean_no_input = Ad_concat*x_0;
x_mean_our_method = x_mean_no_input;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%
iter = 1;
start_time = tic;
while iter <= iter_max 

    % update collision avoid gradient
    [norm_approx_ab, norm_approx_gradient_ab] = ...
            update_g(x_mean_our_method(:,1), x_mean_our_method(:,2), P, Bd_concat, means_ab, time_horizon);
    [norm_approx_ac, norm_approx_gradient_ac] = ...
            update_g(x_mean_our_method(:,1), x_mean_our_method(:,3), P, Bd_concat, means_ac, time_horizon);
    [norm_approx_bc, norm_approx_gradient_bc] = ...
            update_g(x_mean_our_method(:,2), x_mean_our_method(:,3), P, Bd_concat, means_bc, time_horizon);
    
    cvx_begin quiet
        variable U(3 * time_horizon, 3);
        variable x_mean_our_method(6 * time_horizon, 3);

        % 2 vehicle collision avoid
        variable slack(time_horizon, 3) nonnegative;
        
        % targe set
        variable lambda(n_lin_state, 3) nonnegative;
        variable lambda_temp(n_lin_state, 3) nonnegative;
        
        % cost params
        variable sum_slack(1,1);
        variable quad_input_cost(1,1);

        minimize (tau(iter)*sum_slack + quad_input_cost)
        subject to
            %----------------------------
            % cost variables
            %----------------------------
            sum_slack == sum(vec(slack)); 
                      
            quad_input_cost >=  vec(U)'*vec(U);
            
            %----------------------------
            % linear equations defining the state
            %----------------------------
            x_mean_our_method == x_mean_no_input + Bd_concat * U;

            %----------------------------
            % u \in \mathcal{U} 
            %----------------------------
            for i = 1:3
                input_space_A * U(:,i) <= input_space_b;
            end
            
            %----------------------------
            % colission avoidance constraint (intervehicle)
            %----------------------------
%             vec(slack) >= 0;
            for k = 1:time_horizon
                hat_lam * norm( vars_ab((k-1)*4+(1:4),:) * [(x_mean_our_method(6*(k-1)+[1:3], 1) - x_mean_our_method(6*(k-1)+[1:3],2)); 1] ) - ...
                    norm_approx_ab(k) - norm_approx_gradient_ab(k,:) * [U(:,1) - U_p(:,1);U(:,2) - U_p(:,2)] - ...
                    slack(k,1) + r^2 <= 0;
                hat_lam * norm( vars_ac((k-1)*4+(1:4),:) * [(x_mean_our_method(6*(k-1)+[1:3], 1) - x_mean_our_method(6*(k-1)+[1:3],3)); 1] ) - ...
                    norm_approx_ac(k) - norm_approx_gradient_ac(k,:) * [U(:,1) - U_p(:,1);U(:,3) - U_p(:,3)] - ...
                    slack(k,1) + r^2 <= 0;
                hat_lam * norm( vars_bc((k-1)*4+(1:4),:) * [(x_mean_our_method(6*(k-1)+[1:3], 2) - x_mean_our_method(6*(k-1)+[1:3],3)); 1] ) - ...
                    norm_approx_bc(k) - norm_approx_gradient_bc(k,:) * [U(:,2) - U_p(:,2);U(:,3) - U_p(:,3)] - ...
                    slack(k,1) + r^2 <= 0;
            end
            
            %----------------------------
            % terminal state constraint
            %----------------------------
            
            vec(lambda) >= lam_min;
             
            % mean in shrunk target set
            target_set_all_A * (x_mean_our_method(end-5:end, 1) + Wd_concat(end-5:end,:)*W_a_mean) + ...
                 sqrt(diag(target_set_all_A * Wd_concat(end-5:end,:) * W_a_var * Wd_concat(end-5:end,:)' * target_set_all_A')) .* lambda(:,1) - target_set_all_b(:,1) <= 0;
            target_set_all_A * (x_mean_our_method(end-5:end, 2) + Wd_concat(end-5:end,:)*W_b_mean) + ...
                 sqrt(diag(target_set_all_A * Wd_concat(end-5:end,:) * W_b_var * Wd_concat(end-5:end,:)' * target_set_all_A')) .* lambda(:,2) - target_set_all_b(:,2) <= 0;
            target_set_all_A * (x_mean_our_method(end-5:end, 3) + Wd_concat(end-5:end,:)*W_c_mean) + ...
                 sqrt(diag(target_set_all_A * Wd_concat(end-5:end,:) * W_c_var * Wd_concat(end-5:end,:)' * target_set_all_A')) .* lambda(:,3) - target_set_all_b(:,3) <= 0;
            
            for i = 1:(n_lin_state)
                for j = 1:3
                    lambda_temp(i,j) >= pow_func_m .* lambda(i,j) + pow_func_c;
                end
            end
            sum(vec(lambda_temp)) <= 1-safety_target;
            
    cvx_end
    
    % update Costs
    input_cost_our_method(iter+1) = quad_input_cost;
    lambda_sum_our_method(iter+1) = sum_slack;
    total_cost_our_method(iter+1) = cvx_optval;

    % calculate convergence criteria
    conv_check = abs(input_cost_our_method(iter+1) - input_cost_our_method(iter) + tau(iter)*(lambda_sum_our_method(iter+1) - lambda_sum_our_method(iter)));

    % print statistics
    if ~ quiet
        fprintf('iteration: %d ', iter);
        fprintf('\t %f', cvx_optval);
        fprintf('\t %e', conv_check); 
        fprintf('\t %e', lambda_sum_our_method(iter+1));
        fprintf('\t %f', toc(start_time));
        fprintf('\t %s \n', cvx_status);
    end

    % check for solved status
    if strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved')
        % check for convergence
        if (conv_check <= epsilon_dc) && (lambda_sum_our_method(iter+1) <= epsilon_lambda)                 
           break
        end

        % if not converged update previous answer to current answer
        U_p = U;

    % if results are NaN break before error
    elseif strcmpi(cvx_status, 'Failed') || strcmpi(cvx_status, 'Infeasible')
        break
    end

    % update itteration number
     iter = iter + 1;
end
total_time_our_method = toc(start_time);

% make k not less than or equal to iter_max
iter = min(iter, iter_max);


%%%%%%%%%%%%%%%%%%%%%%%%%%
% print some useful information
%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n%s in ', cvx_status);
fprintf('%i itterations \n', iter);
fprintf('Computation time (sec): %f \n', total_time_our_method);
fprintf('Total Cost: %f \n', total_cost_our_method(iter+1));
fprintf('Slack Cost: %f \n', lambda_sum_our_method(iter+1));
fprintf('Input Cost: %f \n', input_cost_our_method(iter+1));

if strcmpi(cvx_status, 'Failed') || strcmpi(cvx_status, 'Infeasible')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify probabilities
%%%%%%%%%%%%%%%%%%%%%%%%%%

[P_target_our_method, P_col_our_method] = verify(x_mean_our_method, Wd_concat, time_horizon, target_sets, r, 10000)
