function [W_mean, W_var] = get_sample_stats_1v(W)
    c = size(W,2);
    W_mean = mean(W,2);
    W_var = zeros(size(W,1));
    for i = 1:c
        W_i = W(:,i) - W_mean;
        W_var = W_var + W_i*W_i';
    end
    W_var = W_var ./c;
end
