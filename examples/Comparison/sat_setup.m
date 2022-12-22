chief_semi_major   = 42164140; % m
chief_eccentricity  = 0;        % unitless 
chief_inclination   = pi / 18;  % rad (10*)
chief_argument_peri = 0;        % rad
chief_raan          = 0;        % rad
chief_true_anom     = pi / 2;   % rad (90*)

dep_1_semi_major    = chief_semi_major + 80;        % m
dep_1_eccentricity  = 0;                            % unitless 
dep_1_inclination   = pi / 18 * (1 - 2e-6);         % rad (10*)
dep_1_argument_peri = 0;                            % rad
dep_1_raan          = 0;                            % rad
dep_1_true_anom     = chief_true_anom * (1 + 2e-7); % rad 

dep_2_semi_major    = chief_semi_major + 85;        % m
dep_2_eccentricity  = 0;                            % unitless 
dep_2_inclination   = pi / 18;                      % rad (10*)
dep_2_argument_peri = 0;                            % rad
dep_2_raan          = 0;                            % rad
dep_2_true_anom     = chief_true_anom * (1 - 1e-7); % rad 

dep_3_semi_major    = chief_semi_major + 87;        % m
dep_3_eccentricity  = 0;                            % unitless 
dep_3_inclination   = pi / 18 * (1 + 1e-6);         % rad (10*)
dep_3_argument_peri = 0;                            % rad
dep_3_raan          = 0;                            % rad
dep_3_true_anom     = chief_true_anom * (1 + 1e-7); % rad 


%%%%%%%%%%%%%%%%
% orbital elements to ECI
%%%%%%%%%%%%%%%%

x_0_chief = kepler2eci(chief_semi_major, chief_eccentricity, chief_inclination, chief_argument_peri, chief_raan, chief_true_anom);
x_0_dep_1 = kepler2eci(dep_1_semi_major, dep_1_eccentricity, dep_1_inclination, dep_1_argument_peri, dep_1_raan, dep_1_true_anom);
x_0_dep_2 = kepler2eci(dep_2_semi_major, dep_2_eccentricity, dep_2_inclination, dep_2_argument_peri, dep_2_raan, dep_2_true_anom);
x_0_dep_3 = kepler2eci(dep_3_semi_major, dep_3_eccentricity, dep_3_inclination, dep_3_argument_peri, dep_3_raan, dep_3_true_anom);

clearvars chief_* dep_*