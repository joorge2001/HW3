function F = srp(t)
% -------------------------------------------------------------------------
% Function to compute the force F, on the sail due to the
% solar radiation pressure when it at a time is t since the solstice 
% ignoring eclipses.

% Inputs: 
%     - t: time elapsed since summer solstice [s]

% Outputs:
%     - F: force produced on the sail due to the solar radiation 
%          pressure [kN]
% -------------------------------------------------------------------------

% Data
rE = 6371; % [km]
%COE AT SUMMER SOLSTICE IN EQUATORIAL ECI RF
zeta_p = 1800; % [km] altitude at pericenter
e = 0.1; % eccentricity 
a = (zeta_p+rE)/(1-e); % [km] semi-major axis
raan = deg2rad(295); % [rad] 
i = deg2rad(33); % [rad]
omega = deg2rad(198); % [rad] 
theta = deg2rad(48); % [rad]

COE = [a,e,i,raan,omega,theta];


P = 4.560e-6; % [N/m^2]
A = 20; % [m^2]
R = 1; % Maximal reflection coefficient

muE = 3.986e5; % [km3/s2]
muS = 1.32712440018e11; % [km3/s2]

a_e = 149597870.7; % [km] we use the value of 1au to set the semi-major axis of the earth's
% orbit

tau_e = 2*pi*sqrt(a_e^3/muS); % [s] period of the Earth's orbit around the Sun
tau = 2*pi*sqrt(COE(1)^3/muE); % [s] period of the satellite's orbit around the Earth

n = 2*pi / tau; % [rad/s] mean motion of the satellite around the earth
n_e = 2*pi / tau_e; % [rad/s] mean motion of the earth around the sun

[r0, v0] = COE2rv(COE(1), COE(2), COE(3), COE(4), COE(5), COE(6) + n*t,muE); % COE to position and velocity vector;

obliquity = deg2rad(23.44);

[r1_ec,~] = EQ2EC(r0,v0,obliquity);


% Calculate the position of the Earth at time t since the solstice
r_earth = [-sin(n_e*t), cos(n_e*t), 0] * a_e; % [km] vector from sun to earth
% x term is negative due to counterclockwise orbit

r_sat_sun = r1_ec' + r_earth;

% Calculate the outward pointing normal to the sail
n_v = -r_sat_sun / norm(r_sat_sun); % The sail is always pointed towards the Sun so n_v = u_v

% Calculate the angle between u_v and n_v
theta = acos(dot(n_v, n_v));

% Calculate force F due to solar radiation pressure
F = -P*A*cos(theta)*((1-R)*n_v + 2*R*cos(theta)*n_v);

% Convert to kN
F = F * 1e-3;

end


