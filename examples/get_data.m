function [data] = get_data(filename, n)
    if isfile(filename)
        data = readmatrix(filename);
        return        
    end
    data = zeros(48, n);
    cov = diag([10^-6, 10^-6, 10^-6, 5*10^-8, 5*10^-8, 5*10^-8]);
    for i = 1:n
        data(:,i) = mvtrnd(kron(eye(8), cov),1)'.*10^-6 + exprnd(ones(48,1)*10^-4);
    end
    writematrix(data, filename);
end