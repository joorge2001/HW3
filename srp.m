function F = srp(t)
    % -------------------------------------------------------------------------
    % Function to compute the force F, on the sail due to the
    % solar radiation pressure when it at a time is t since the solstice 
    % ignoring eclipses.
    
    % Inputs: 
    %     - t: time since solstice [s]

    % Outputs:
    %     - F: force produced on the sail due to the solar radiation 
    %          pressure [N]
    % -------------------------------------------------------------------------

% Data
h = 6.62607015e-34; % [J/Hz] Planck's constant
S_e = 1367; % [W/m^2] goes with 1/r^2
P = 4.560e-6; % [N/m^2]
A = 20; % [m^2]
R = 1; % Maximal reflection coefficient
u_v = [] %unit vector that points to the sun (measured from earth as it is considered
% to be irrelevant if the earth or the spacecraft is the origin of the
% vector
n_v = [] %outward pointing normal (to sail) unit vector n that forms an 
% angle theta with u_v
% The orbit is circular: u_v = - n_v and theta = 0;

F = -P*A*cos(theta)*((1-R)*u_v + 2*R*cos(theta)*n_v);


