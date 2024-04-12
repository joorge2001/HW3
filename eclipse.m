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

% Data
rE = 6371; % [km] Earth's radius
a_e = 149597870.7;
muS = 1.32712440018e11; % [km3/s2]
tau_e = 2*pi*sqrt(a_e^3/muS); % [s] period of the Earth's orbit

% Calculate the position vector of the earth and satellite seen in the 
% sun centric reference frame
r_earth = [-sin(2*pi*t/tau_e), cos(2*pi*t/tau_e), 0] * a_e; % [km] vector from sun to earth

r_sat_sun = r' + r_earth;

% Calculate the angle between the position vectors of the Earth and the satellite relative to the Sun
angle = abs(acos(dot(r_sat_sun, r_earth) / (norm(r_sat_sun) * norm(r_earth))));

% Calculate the critical angle
critical_angle = abs(atan(rE/a_e));

disp(angle);

d = (norm(r_sat_sun)*cos(angle));

disp(d);
if (angle < critical_angle) && (d > a_e)
    flag = true;
else
    flag = false;
end

end