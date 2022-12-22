data_a = readmatrix('./Data Collection/a_dist_test.dat');
data_b = readmatrix('./Data Collection/b_dist_test.dat');
data_c = readmatrix('./Data Collection/c_dist_test.dat');

samples = size(data_a,2);

target_checker = ones(samples, 1);
col_avoid_chief = ones(samples, 1);
col_avoid_dep = ones(samples, 1);

for i = 1 :samples
    dep1 = Ex_A_no + C_concat * Ua + D_concat * data_a(:,i);
    dep2 = Ex_B_no + C_concat * Ub + D_concat * data_b(:,i);
    dep3 = Ex_C_no + C_concat * Uc + D_concat * data_c(:,i);
    
    % check terminal sets
    if any(target_set_a.A * dep1(end-5:end) > target_set_a.b)
        target_checker(i) = 0;
    end
    if any(target_set_b.A * dep2(end-5:end) > target_set_b.b)
        target_checker(i) = 0;
    end
    if any(target_set_c.A * dep3(end-5:end) > target_set_c.b)
        target_checker(i) = 0;
    end
    
    for j = 1:time_steps
        % check collision avoid betwen deputies
        if (norm(Sk(:,:,j) * (dep1 - dep2)) < r)
            col_avoid_dep(i) = 0;
        end
        if (norm(Sk(:,:,j) * (dep1 - dep3)) < r)
            col_avoid_dep(i) = 0;
        end
        if (norm(Sk(:,:,j) * (dep2 - dep3)) < r)
            col_avoid_dep(i) = 0;
        end
           
        % check collision avoid with chief
        if (norm(Sk(:,:,j) * dep1) < r)
            col_avoid_chief(i) = 0;
        end
        if (norm(Sk(:,:,j) * dep2) < r)
            col_avoid_chief(i) = 0;
        end
        if (norm(Sk(:,:,j) * dep3) < r)
            col_avoid_chief(i) = 0;
        end
    end
    
end

p_target = mean(target_checker)
p_chief = mean(col_avoid_chief)
p_dep = mean(col_avoid_dep)