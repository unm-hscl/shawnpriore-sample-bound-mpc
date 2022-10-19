function [means, vars] = get_sample_stats_2v(W, P)
    means = zeros(4*8,5);
    vars = zeros(4*8,4);
    for i = 1:8
        index = (i-1)*6 + [1:6];
        W_k = P*W(index,:);
        m_hat = mean(W_k,2);
        m_tilde = 0;
        for j = 1:100
           m_tilde = m_tilde + W_k(:,j)'*W_k(:,j)/100;
        end
        s_hat = zeros(3);
        s_tilde = 0;
        s_breve = zeros(3,1);
        for j = 1:100
            s_hat = s_hat + (W_k(:,j) - m_hat) * (W_k(:,j) - m_hat)' ./100;
            s_tilde = s_tilde + (W_k(:,j)'*W_k(:,j) - m_tilde)^2 /100;
            s_breve = s_breve + (W_k(:,j) - m_hat) .* (W_k(:,j)'*W_k(:,j) - m_tilde) ./100;
        end
        means((i-1)*4+[1:4],1:4) = chol([eye(3), m_hat; m_hat', m_tilde]);
        means((i-1)*4+[1:3],5) = m_hat;
        vars((i-1)*4+[1:4],:) = chol([4*s_hat, 2*s_breve; 2*s_breve', s_tilde]);
    end
end