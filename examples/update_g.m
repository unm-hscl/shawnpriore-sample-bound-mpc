function [norm_est, gradient_est] = update_g(x_1, x_2, S, Bd, means, time_horizon)
    norm_est = zeros(time_horizon, 1);
    grad_temp = zeros(time_horizon, size(Bd,2));

    mu = x_1 - x_2;
    for k = 1:time_horizon
        index = 6*(k-1) + (1:6);
        y_k = S * mu(index);
        norm_est(k) = norm(means(4*(k-1)+(1:4),1:4)*[y_k;1])^2;
                
        % get indexed rows of controlability matrix
        Cu_i = Bd(index, :);
        
        % calculate gradient of norm
        grad_temp(k,:) =  2 * y_k' * S * Cu_i + 2 * means((k-1)*4+[1:3],5)' * S * Cu_i;
    end
    gradient_est = [grad_temp, -grad_temp];   
end