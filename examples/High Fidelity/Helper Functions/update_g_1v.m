function [norm_est, gradient_est] = update_g_1v(x_1, S, Bd, W_means_mat, time_horizon)
    norm_est = zeros(time_horizon, 1);
    gradient_est = zeros(time_horizon, size(Bd,2));
    
    state_size = size(S,2);
    
    Sk = zeros(3, state_size * time_horizon, time_horizon);
    for i = 1:time_horizon
        Sk(:, state_size*(i-1) + (1:state_size), i) = S;
    end

    for k = 1:time_horizon
        x_k = Sk(:,:,k) * x_1;
        norm_est(k) = [x_k; 1]' * W_means_mat(:,:,k) * [x_k; 1];
        gradient_est(k,:) =  2 * [x_k; 1]' * W_means_mat(:,1:3,k) * Sk(:,:,k) * Bd;
    end
end