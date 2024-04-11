function flag = eclipse(r,t)
% -------------------------------------------------------------------------
% Function to check if the position vector is inside the eclipse of the
% Earth at a certain time since the solstice.

% Inputs: 
%     - r: position vector of the satellite in ECI-ecliptic reference
%          frame [km] [3x1]
%     - t: time elapsed since summer solstice [s]

% Outputs:
%     - flag: boolean value indicating if the satellite is in the eclipse
%     of the Earth.
%       - flag = true when it is inside
%       - flag = false otherwise
% -------------------------------------------------------------------------

a_e = 149597870.7;
muS = 1.32712440018e11; % [km3/s2]
tau_e = 2*pi*sqrt(a_e^3/muS); % [s] period of the Earth's orbit

r_earth = [cos(2*pi*t/tau_e), sin(2*pi*t/tau_e), 0] * a_e; % [km] vector from sun to earth

r_sat_sun = r0_ec' + r_earth;

if norm(r_sat_sun) < a_e 
    flag = true;
else
    flag = false;
end



end