function [means, vars] = get_sample_stats_2v(S, Dk_concat, W1, W2, time_steps)
    samples = size(W1,2);
    means = zeros(size(S,1)+1, size(S,1)+1, time_steps);
    vars = zeros(size(S,1)+1, size(S,1)+1, time_steps);
    for i = 1:time_steps
        index = (i-1)*6 + [1:6];
        Dk = Dk_concat(index, :);
        
        Zk = S * Dk * (W1 - W2);
        
        EX_Zk = mean(Zk,2);
        EX_ZkTZk = 1/samples * trace(Zk'*Zk);
        
        means(:,:,i) = [eye(size(S,1)), EX_Zk; EX_Zk', EX_ZkTZk];
        
        VAR_Zk = zeros(size(S,1));
        VAR_ZkTZk = 0;
        COV_Zk_ZkTZk = zeros(size(S,1),1); 
        
        for sample = 1:samples
            Zk_i = Zk(:,sample) - EX_Zk;
            Zk_iTZk_i = Zk_i'*Zk_i - EX_ZkTZk;
            
            VAR_Zk = VAR_Zk + (Zk_i*Zk_i' ./ samples);
            VAR_ZkTZk = VAR_ZkTZk + (Zk_iTZk_i^2 / samples);
            COV_Zk_ZkTZk = COV_Zk_ZkTZk + (Zk_i .* Zk_iTZk_i ./ samples);
        end
        vars(:,:,i) = chol([4*VAR_Zk, 2*COV_Zk_ZkTZk; 2*COV_Zk_ZkTZk', VAR_ZkTZk]);
        
    end
end