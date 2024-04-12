function [COE] = hw3data()

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
end