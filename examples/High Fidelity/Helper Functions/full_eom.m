function xdot = full_eom(t, x, JD)
    xdot = zeros(6,1);

    % constants
    mu      = 3.986004418e14;   % m^3/sec^2
    mu_sun  = 1.32712440042e20; % m^3/sec^2
    mu_moon = 4.9048695e12;     % m^3/sec^2
    rad_earth = 6.3781e6; % m
    J2_coef = 1.082e-3;
    deg =  pi / 180;
    per_asec = 1/3600;
    epsilon = 23.43929111 * deg;

    % split vec
    x_sat = x(1:3);
    xdot_sat = x(4:6);

    %compute distance to center of earth
    x_sat_norm = (x_sat'*x_sat)^(0.5);

    % GM accel
    simple_accel    = - x_sat .* (mu /  x_sat_norm^3);

    % J2 harmonic accel
    J2_accel        = ((3 * mu * J2_coef * rad_earth^2) / (2 * x_sat_norm^5)) * ((5* x(3)^2 / x_sat_norm^2)*ones(3,1) - [1;1;3]) ./ x_sat;

    % sun accel
    T = (JD - 2451545) / 36525;
    M = (357.5256 + 35999.049 * T) * deg;

    sun_long = 282.94 * deg + M + (6892 * sin(M) + 72 * sin(2*M)) * deg * per_asec;
    sun_dist = (149.619 - 2.499 * cos(M) - 0.021 * cos(2*M) ) * 1e9; % m
    x_sun = sun_dist * [ cos(sun_long); sin(sun_long) * cosd(epsilon); sin(sun_long) * sind(epsilon)];

    sun_diff = x_sun - x_sat;
    sun_diff_norm = (sun_diff'*sun_diff)^(0.5);

    sun_accel = mu_sun*(sun_diff/sun_diff_norm^3 - x_sun/sun_dist^3);

    % moon accel
    L0  = (218.31617 + 481267.88088 * T - 1.3972 * T) * deg;
    l   = (134.96292 + 477198.86753 * T) * deg;
    lp  = (357.52543 +  35999.04944 * T) * deg;
    F   = ( 93.27283 + 483202.01873 * T) * deg;
    D   = (297.85027 + 445267.11135 * T) * deg;

    moon_long = L0 ...
        + deg * per_asec * (22640 * sin(l) + 769 * sin(2 * l) ...
            - 4586 * sin(l - 2 * D) + 2370 * sin(2 * D) ...
            - 668 * sin(lp) - 413 * sin(2 * F) ...
            - 212 * sin(2 * l - 2 * D) - 206 * sin(l + lp - 2 * D) ...
            + 192 * sin(l + 2 * D) - 165 * sin(lp - 2 * D) ...
            + 148 * sin(l - lp) - 125 * sin(D)...
            - 110 * sin(l + lp) - 55 * sin(2 * F - 2 * D) );
    moon_lat = deg * per_asec ...
        * (18520 * sin(F + moon_long - L0 + (413 * sin(2 * F) + 541 * sin(lp)) * deg * per_asec) ...
            - 526 * sin(F - 2*D) + 44 * sin(l + F - 2*D) ...
            - 31 * sin(-l + F - 2 * D) - 25 * sin(-2 * l + F) ...
            - 23 * sin(lp + F - 2 * D) + 21 * sin(-l + F) ...
            + 11 * sin(-lp + F - 2 * D));
    moon_dist = (385000 - 20905 * cos(l) - 3699 * cos(2 * D - l) ...
        - 2956 * cos(2 * D) - 570 * cos(2 * l) + 246 * cos(2 * l - 2 * D) ...
        - 205 * cos(lp - 2 * D) - 171 * cos(l + 2 * D) - 152 * cos(l + lp - 2 * D)) * 1e3; % m

    ecliptic_rotation = [1, 0, 0; 0, cos(-epsilon), sin(-epsilon); 0, -sin(-epsilon), cos(-epsilon)];
    x_moon = moon_dist * ecliptic_rotation * [ cos(moon_long)*cos(moon_lat); sin(moon_long)*cos(moon_lat); sin(moon_lat)];

    moon_diff = x_moon - x_sat;
    moon_diff_norm = (moon_diff'*moon_diff)^(0.5);

    moon_accel = mu_moon * (moon_diff/moon_diff_norm^3 - x_moon/moon_dist^3);

    % eom
    xdot(1:3) = xdot_sat;
    xdot(4:6) = simple_accel + J2_accel + sun_accel + moon_accel;
   
end