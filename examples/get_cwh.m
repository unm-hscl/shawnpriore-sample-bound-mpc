function [Ad_concat, Bd_concat, Wd_concat] = get_cwh(time_horizon)

    sampling_period              = 30;                                              % sec
    orbital_radius               = (2000 + 6378.1) * 1000;                          % m
    gravitational_constant       = 6.673e-11;                                       % m^3 kg^-1 sec^-2
    celestial_mass               = 5.9472e24;                                       % kg
    gravitational_body           = gravitational_constant * celestial_mass;         % m^3 sec^-2
    orbit_ang_vel                = sqrt(gravitational_body / orbital_radius^3);     % rad sec^-2

    % Continuous-time LTI CWH unforced dynamics e^{At}
    e_power_At = @(t) [ 
        4 - 3 * cos(orbit_ang_vel * t), 0, 0, (1/orbit_ang_vel) * sin(orbit_ang_vel * t), (2/orbit_ang_vel) * (1 - cos(orbit_ang_vel * t)), 0; 
        6 * (sin(orbit_ang_vel * t) - orbit_ang_vel * t), 1, 0, -(2/orbit_ang_vel) * (1 - cos(orbit_ang_vel * t)), (1/orbit_ang_vel) * (4*sin(orbit_ang_vel * t) - 3*orbit_ang_vel * t), 0; 
        0, 0, cos(orbit_ang_vel * t), 0, 0, (1/orbit_ang_vel) * sin(orbit_ang_vel * t); 
        3 * orbit_ang_vel * sin(orbit_ang_vel * t), 0, 0, cos(orbit_ang_vel * t), 2 * sin(orbit_ang_vel * t), 0; 
        -6 * orbit_ang_vel * (1 - cos(orbit_ang_vel * t)), 0, 0, -2 * sin(orbit_ang_vel * t), 4 * cos(orbit_ang_vel * t) - 3, 0;
        0, 0, -orbit_ang_vel * sin(orbit_ang_vel * t), 0, 0, cos(orbit_ang_vel * t);
        ];

    % Discrete-time system is Phi(T_s) for sampling time T_s since the system is time-invariant
    Ad = e_power_At(sampling_period);
    
    % Impulse control
    Bd = Ad*[zeros(3); eye(3)];
    
    % Concat matrix maker
    Ad_concat = zeros(size(Ad, 1)*time_horizon, size(Ad, 2));
    Bd_concat = zeros(size(Bd, 1)*time_horizon, size(Bd, 2)*time_horizon);
    Wd_concat = zeros(size(Ad, 1)*time_horizon);
    for i = 0:(time_horizon-1)
        Ad_concat(size(Ad, 1)*i + [1:size(Ad, 1)], :) = Ad^(i+1);
    end
    for i = 0:(time_horizon-1)
        for j = 0:i
            Bd_concat(size(Bd, 1)*i + [1:size(Bd, 1)], size(Bd, 2)*j + [1:size(Bd, 2)]) = Ad^(i-j) * Bd;
            Wd_concat(size(Ad, 1)*i + [1:size(Ad, 1)], size(Ad, 2)*j + [1:size(Ad, 2)]) = Ad^(i-j);
        end
    end
end