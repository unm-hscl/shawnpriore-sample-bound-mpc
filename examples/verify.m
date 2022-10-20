function [P_target, P_collision] = verify(mean_sat, Wd_concat, time_horizon, target_sets, r, samples)
%   VERIFY Summary of this function goes here
%   Detailed explanation goes here

    vehicles = size(mean_sat, 2);
    combinations = vehicles * (vehicles-1) / 2;
        
    disturbed_mean_sats = zeros([size(mean_sat),samples]);
    cov = diag([10^-5, 10^-5, 10^-5, 5*10^-8, 5*10^-8, 5*10^-8]);
    
    rng('default');
    for i = 1:vehicles
        for j = 1:samples
            disturbed_mean_sats(:,i,j) =  mvtrnd(kron(eye(8), cov),1)'.*10^-5 + exprnd(ones(48,1)*10^-4);
        end
    end
    
    for i = 1:samples
        disturbed_mean_sats(:,:,i) = mean_sat +  Wd_concat * disturbed_mean_sats(:,:,i);
    end
   
    in_target = zeros(vehicles, samples);
    for i = 1:vehicles
        for j = 1:samples
            in_target(i, j) = target_sets(i).contains( disturbed_mean_sats(end-5:end, i, j) );
        end
    end
    P_target = sum(sum(in_target, 1) == vehicles) / samples;
    
    
    collision_sats = zeros(combinations, time_horizon, samples);
    for i = 1:(vehicles-1)
        for j = (i+1):vehicles
            index = (i-1)*(vehicles-1-i/2) + j-1;
            for t = 1:time_horizon
                time_index = 6*(t-1) + [1:3];
                for k = 1:samples
                    collision_sats(index, t, k) = ( norm( disturbed_mean_sats(time_index, i, k) - disturbed_mean_sats(time_index, j, k)) >= r);
                end
            end
        end
    end
    P_collision = sum(sum(sum(collision_sats, 1) == combinations, 2) == time_horizon) / samples;
 end

