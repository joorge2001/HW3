function [r0_ec,v0_ec] = EQ2EC(r0,v0,obliquity)
    % -------------------------------------------------------------------------
    % Function to compute position and velocity vectors in ECI-ecliptic, from 
    % a position and velocity vector in ECI-equatorial.
    
    % Inputs: 
    %     - r0: position vector in ECI-eq reference frame [3x1]
    %     - v0: velocity vector in ECI-eq reference frame [3x1]
    %     - obliquity = angle between the bodie's equatorial plane and its 
    %                   orbital plane.
    %     
    % Outputs:
    %     - r0_ec: position vector in ECI-ec reference frame [3x1]
    %     - v0_ec: velocity vector in ECI-ec reference frame [3x1]
    % -------------------------------------------------------------------------

% Rotation matrix
R = [1, 0, 0;
     0, cos(obliquity), sin(obliquity);
     0, -sin(obliquity), cos(obliquity)];

% Convert position and velocity to ecliptic ECI reference frame
r0_ec = R * r0;
v0_ec = R * v0;