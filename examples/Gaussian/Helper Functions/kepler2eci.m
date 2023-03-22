function vec = kepler2eci(axi, ecc, inc, arg, raan, tru)

mu = 3.986004418e14; % m^3/sec^2

p = axi*(1-ecc ^2);
r_0 = p / (1 + ecc * cos(tru));

x = r_0 * cos(tru);
y = r_0 * sin(tru);

Vx_ = -(mu/p)^(1/2) * sin(tru);
Vy_ = (mu/p)^(1/2) * (ecc + cos(tru));

X = (cos(raan) * cos(arg) - sin(raan) * sin(arg) * cos(inc)) * x + (-cos(raan) * sin(arg) - sin(raan) * cos(arg) * cos(inc)) * y;
Y = (sin(raan) * cos(arg) + cos(raan) * sin(arg) * cos(inc)) * x + (-sin(raan) * sin(arg) + cos(raan) * cos(arg) * cos(inc)) * y;
Z = (sin(arg) * sin(inc)) * x + (cos(arg) * sin(inc)) * y;
Vx = (cos(raan) * cos(arg) - sin(raan) * sin(arg) * cos(inc)) * Vx_ + (-cos(raan) * sin(arg) - sin(raan) * cos(arg) * cos(inc)) * Vy_;
Vy = (sin(raan) * cos(arg) + cos(raan) * sin(arg) * cos(inc)) * Vx_ + (-sin(raan) * sin(arg) + cos(raan) * cos(arg) * cos(inc)) * Vy_;
Vz = (sin(arg) * sin(inc)) * Vx_ + (cos(arg) * sin(inc)) * Vy_;

vec = [X; Y; Z; Vx; Vy; Vz];
end
