function x_cart = polar2cart(x_polar,theta)
% -------------------------------------------------------------------------
% Function to compute position and velocity vectors in ECI, from a
% position and velocity vector in SZI.

% Inputs: 
%     - x_polar: vector in polar basis [3x1]
%     - theta: true anomaly [rad]
%     
% Outputs:
%     - x_cart: vector in cartesian basis [3x1]
% -------------------------------------------------------------------------

% Rotation matrix
R = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];

% Final vector
x_cart = R*x_polar;

end