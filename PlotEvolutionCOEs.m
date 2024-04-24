function [] = PlotEvolutionCOEs(t_span, X)
% -------------------------------------------------------------------------
% Function to go from classical orbital elements to position and velocity
% vectors in ECI reference frame.
%
% Inputs: 
%     - t_span: time span obtained from ode45
%     - X: state vector obtained from ode45
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

muE = 3.986e5; % [km3/s2]

% Initialize Orbital Elements
a = 0*t_span;
e = 0*t_span;
i = 0*t_span;
RAAN = 0*t_span;
omega = 0*t_span;
theta = 0*t_span;

for k = 1:length(t_span)
    rr = X(k, 1:3);
    vv = X(k, 4:6);
    
    % Obtain orbital elements
    [a(k),e(k),i(k),RAAN(k),omega(k),theta(k)] = rv2COE(muE,rr,vv);
    
    i(k) = rad2deg(i(k));
    RAAN(k) = rad2deg(RAAN(k));
    omega(k) = rad2deg(omega(k));
    theta(k) = rad2deg(theta(k));
    
end

% Plot the evolution of the Orbital Elements over time
figure
plot(t_span, a);
title('Semimajor Axis: $a$','Interpreter', 'latex')
xlabel('time [s]','Interpreter', 'latex')
ylabel('a [km]','Interpreter', 'latex')

figure
plot(t_span, e);
title('Eccentricity: $e$','Interpreter', 'latex')
xlabel('time [s]','Interpreter', 'latex')
ylabel('e','Interpreter', 'latex')

figure
plot(t_span, i);
title('Inclination: $i$','Interpreter', 'latex')
xlabel('time [s]','Interpreter', 'latex')
ylabel('$i [^\circ]$','Interpreter', 'latex')

figure
plot(t_span, RAAN);
title('RAAN: $\Omega$','Interpreter', 'latex')
xlabel('time [s]','Interpreter', 'latex')
ylabel('$\Omega [^\circ]$','Interpreter', 'latex')

figure
plot(t_span, omega);
title('Argument of Periapsis: $\omega$','Interpreter', 'latex')
xlabel('time [s]','Interpreter', 'latex')
ylabel('$\omega [^\circ]$','Interpreter', 'latex')

figure
plot(t_span, theta);
title('True Anomaly: $\theta$','Interpreter', 'latex')
xlabel('time [s]','Interpreter', 'latex')
ylabel('$\theta [^\circ]$','Interpreter', 'latex')


end