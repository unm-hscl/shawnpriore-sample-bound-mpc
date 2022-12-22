function xdot = simple_eom(t,x)
    xdot = zeros(6,1);

    % constants
    mu = 3.986e14; %  [m^3/sec^2]

    % calculations
    x_t = x(1:3);
    xdot_t = x(4:6);

    x_t_norm = (x_t'*x_t)^(0.5);

    % eom
    xdot(1:3) = xdot_t;
    xdot(4:6) = - x_t .* (mu /  x_t_norm^3);

end