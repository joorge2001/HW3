function x_ECI = Perifocal2ECI(x_p,raan,i,omega)
% -------------------------------------------------------------------------
% Function go from Perifocal to ECI reference frame.

% Inputs: 
%     - x_p: vector in perifocal basis [3x1]
%     - raan: right ascension of the ascending node [rad]
%     - i: inclination [rad]
%     - omega: argument of periapsis [rad]
%
% Outputs:
%     - x_ECI: vector in ECI basis [3x1]
% -------------------------------------------------------------------------

% Rotation Matrix

R = [cos(omega)*cos(raan)-cos(i)*sin(omega)*sin(raan), ...
     -sin(omega)*cos(raan)-cos(i)*cos(omega)*sin(raan),...
     sin(i)*sin(raan);...
     cos(omega)*sin(raan)+cos(i)*sin(omega)*cos(raan), ...
     -sin(omega)*sin(raan)+cos(i)*cos(omega)*cos(raan),...
     -sin(i)*cos(raan);...
     sin(i)*sin(omega),...
     sin(i)*cos(omega),...
     cos(i)];
 
x_ECI = R*x_p;
 
end