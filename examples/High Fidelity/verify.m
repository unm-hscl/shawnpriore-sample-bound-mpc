function [p, samples] = verify(traj, D_concat, G, h)
    data = readmatrix('./Data Collection/a_dist_test.dat');
    samples = size(data,2);

    target_checker = zeros(samples, 1);

    for i = 1 :samples
        x = traj + D_concat * data(:,i);
        target_checker(i) = all(G*x <= h);
    end

    p = mean(target_checker);
end