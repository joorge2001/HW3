function [r0,v0] = COE2rv(a,e,i,raan,omega,theta,mu)
% -------------------------------------------------------------------------
% Function to go from classical orbital elements to position and velocity
% vectors in ECI reference frame.
%
% Inputs: 
%     - a: semimajor axis [km]
%     - e: eccentricity
%     - i: inclination [rad]
%     - raan: right ascension of the ascending node [rad]
%     - omega: argument of periapsis [rad]
%     - theta: true anomaly [rad]
%     - mu: gravitational parameter [km3/s2]
%
% Outputs:
%     - r0: position vector in ECI reference frame [3x1]
%     - v0: velocity vector in ECI reference frame [3x1]
% -------------------------------------------------------------------------

p = a*(1-e^2);
h = sqrt(mu*p);

r_norm = p/(1+e*cos(theta)); % trajectory equation
% v_norm = sqrt(2*my/r_norm - mu/a); % vis viva equation

r_p_polar = [r_norm; 0; 0]; % position vector in perifocal ref. frame in polar 
                       % basis

v_r = (mu/h)*e*sin(theta); % radial component of velocity
v_theta = (mu/h)*(1+e*cos(theta)); % theta component of velocity

v_p_polar = [v_r; v_theta; 0]; % velocity vector in perifocal ref. frame 
                               % in polar basis

% Transform vectors from polar basis to cartesian basis
r_p_cart = polar2cart(r_p_polar,theta); % vector still in perifocal frame
v_p_cart = polar2cart(v_p_polar,theta); % vector still in perifocal frame       

% Transform vectors from perifocal to ECI reference frame
r0 = Perifocal2ECI(r_p_cart,raan,i,omega); % vector in ECI frame
v0 = Perifocal2ECI(v_p_cart,raan,i,omega); % vector in ECI frame   

end

