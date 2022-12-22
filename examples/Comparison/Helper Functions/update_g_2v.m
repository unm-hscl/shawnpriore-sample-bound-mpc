function [norm_est, gradient_est] = update_g_2v(x_1, x_2, S, Bd, W_means_mat, time_horizon)
    norm_est = zeros(time_horizon, 1);
    grad_temp = zeros(time_horizon, size(Bd,2));
    
    state_size = size(S,2);
    
    Sk = zeros(3, state_size * time_horizon, time_horizon);
    for i = 1:time_horizon
        Sk(:, state_size*(i-1) + (1:state_size), i) = S;
    end

    mu = x_1 - x_2;
    for k = 1:time_horizon
        x_k = Sk(:,:,k) * mu;
        norm_est(k) = [x_k; 1]' * W_means_mat(:,:,k) * [x_k; 1];
        grad_temp(k,:) =  2 * [x_k; 1]' * W_means_mat(:,1:3,k) * Sk(:,:,k) * Bd;
    end
    gradient_est = [grad_temp, -grad_temp];   
end