function F = srp(t)
% -------------------------------------------------------------------------
% Function to compute the force F, on the sail due to the
% solar radiation pressure when it at a time is t since the solstice 
% ignoring eclipses.

% Inputs: 
%     - t: time elapsed since summer solstice [s]

% Outputs:
%     - F: force produced on the sail due to the solar radiation 
%          pressure [N]
% -------------------------------------------------------------------------

% Data
COE = hw3data();
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


truean = COE(6)+n*t;

while truean > 2*pi
    truean = truean - 2*pi;
end

while truean < 0
    truean = truean + 2*pi;
end

[r0, v0] = COE2rv(COE(1), COE(2), COE(3), COE(4), COE(5), COE(6) + n*t,muE); % COE to position and velocity vector;

obliquity = deg2rad(23.44);

[r0_ec,~] = EQ2EC(r0,v0,obliquity);

flag = eclipse(r0_ec,t); % r0_ec is in the ecliptic frame and points from earth to satellite
disp(flag);


% Calculate the position of the Earth at time t since the solstice
r_earth = [-sin(2*pi*t/tau_e), cos(2*pi*t/tau_e), 0] * a_e; % [km] vector from sun to earth
% x term is negative due to prograde orbit

% Calculate the unit vector pointing to the Sun
u_v = -r_earth / norm(r_earth); % The Sun is in the opposite direction of the Earth

r_sat_sun = r0_ec' + r_earth;

% Calculate the outward pointing normal to the sail
n_v = r_sat_sun / norm(r_sat_sun); % The sail is always pointed towards the Sun

% Calculate the angle between u_v and n_v
theta = acos(dot(u_v, n_v));

% Calculate force F due to solar radiation pressure
F = -P*A*cos(theta)*((1-R)*u_v + 2*R*cos(theta)*n_v);


