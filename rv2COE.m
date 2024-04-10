function [a,e,i,raan,omega,theta] = rv2COE(mu,r0,v0)
% -------------------------------------------------------------------------
% Function to compute classical orbital elements, from a position 
% and velocity vector in a certain reference frame, at the instant of the observation.

% Inputs: 
%     - mu: gravitational parameter [km3/s2]
%     - r0: position vector [3x1]
%     - v0: velocity vector [3x1]
%
% Outputs:
%     - a: semimajor axis [km]
%     - e: eccentricity
%     - i: inclination [rad]
%     - raan: right ascension of the ascending node [rad]
%     - omega: argument of periapsis [rad]
%     - theta: true anomaly [rad]
% -------------------------------------------------------------------------

% To be compliant with the usual range of values used in orbital mechanics,
% conditional statements have been employed to ensure the proper definition
% of the angles in this function.

% Unitary vectors
i_v = [1;0;0];
j_v = [0;1;0];
k_v = [0;0;1];

% Obtain the semimajor axis with the vis-viva equation
r = norm(r0); % [km/s]
v = norm(v0); % [km/s]
a = -mu/2 * 1/((v^2)/2-mu/r); % [km]

% Eccentricity vector and eccentricity
h_v = cross(r0,v0); % [km^2/s]
h = norm(h_v); % [km^2/s]

e_v = (v^2-mu/r).*r0./mu - dot(r0,v0)./mu.*v0;
e = norm(e_v);

% Inclination
cosi = dot(h_v,k_v)/h;
sini = sqrt(1-cosi^2);

if cosi <= 0 || sini <=0
    i = atan2(sini,cosi); % [rad]
else
    i = atan(sini/cosi); % [rad]
end

% Right ascension of the ascending node
n = cross(k_v,h_v)/norm(cross(k_v,h_v));

cosraan = dot(i_v,n);
sinraan = dot(j_v,n);

if cosraan <= 0 || sinraan <=0
    raan = atan2(sinraan,cosraan); % [rad]
else
    raan = atan(sinraan/cosraan); %[rad]
end

if raan <0 
    raan = raan +2*pi;
end

% Argument of periapsis
m = cross(h_v,n)/h;
cosomega = dot(n,e_v)/e; 
sinomega = dot(m,e_v)/e; 

if cosomega <= 0 || sinomega <=0
    omega = atan2(sinomega,cosomega); % [rad]
else
    omega = atan(sinomega/cosomega); % [rad]
end

if omega <0 
    omega = omega +2*pi; 
end

% True anomaly
p = cross(h_v,e_v)/(h*e);
costheta = dot(e_v,r0)/(e*r);
sintheta = dot(p,r0)/(r);

if costheta <= 0 || sintheta <=0
    theta = atan2(sintheta,costheta); % [rad]
else
    theta = atan(sintheta/costheta); % [rad]
end

while theta > 2*pi
    theta = theta - 2*pi;
end

while theta < 0
    theta = theta + 2*pi;
end

end