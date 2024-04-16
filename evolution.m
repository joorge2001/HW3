function dX = evolution(t, X)
% -------------------------------------------------------------------------
% Function calculate the time derivative of the state vector of the
% satellite taking into account the attracton of the Earth and the solar
% radiation pressure (only present outside of the eclipses)
%
% Inputs: 
%     - X: state vector of the satellite in ECI-ecliptic reference
%          frame. It has the form of dX = [rx ry rz vx vy vz], where each
%          of these correspond to the position and velocity components of
%          the satellite. Units are [km] and [km/s]
%     - t: time elapsed since Summer solstice [s]

% Outputs:
%     - dX: time derivative of the state vector with respect to the 
%           ECI-ecliptic reference frame
% -------------------------------------------------------------------------


% Data
m = 100;        % mass of satellite [kg]
muE = 3.986e5;  % [km3/s2]

% Define position and velocity vectors
rr = X(1:3);
vv = X(4:6);
r = norm(rr);

% Compute main acceleration towards the Earth
a_m = -muE/r^3 * rr;

% Compute the force produced by the solar radiation
if eclipse(rr, t) == false
    F_solar = srp(t)';
else
    F_solar = [0; 0; 0];
end


% Compute time derivative of state vector
dX (1:3,1) = vv; 
dX (4:6,1) = a_m + (F_solar/m);

end




