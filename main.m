%%-------- Homework 3: Space Vehicles and Orbital Dynamics --------

% Authors: 
%   - Jaume Manresa Rigo (100461930)
%   - Jorge Cortés de Jesús (100428594)

%%
clear; clc; 
close all;

%% Format settings

set(groot, 'defaultLegendFontSize', 20);
set(groot, 'defaultTextFontSize', 20);
set(groot, 'defaultAxesFontSize', 20);
set(groot, 'defaultAxesLineWidth', 1);
set(groot, 'defaultAxesXMinorTick', 'on');
set(groot, 'defaultAxesYMinorTick', 'on');
set(groot, 'defaultLegendBox', 'off');
set(groot, 'defaultLegendLocation', 'best');
set(groot, 'defaultLineLineWidth', 1);
set(groot, 'defaultLineMarkerSize', 10);
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% 

m = 100; % [kg]
S_ss = 20; % [m^2] always points towards the sun
% orbit around earth

rB = 6371; % [km]
muE = 3.986e5; % [km3/s2]
muS = 1.32712440018e11; % [km3/s2]

%COE AT SUMMER SOLSTICE IN EQUATORIAL ECI RF
zeta_p = 1800; % [km] altitude at pericenter
e = 0.1; % eccentricity 
a = (zeta_p+rB)/(1-e); % [km] semi-major axis
raan = deg2rad(295); % [rad] 
i = deg2rad(33); % [rad]
omega = deg2rad(198); % [rad] 
theta = deg2rad(48); % [rad]

COE = [a,e,i,raan,omega,theta];
name = {'Earth';'Satellite'};
ellipticOrbitPlotter(COE,muE,name,rB,1)

% consider the orbit of the earth around the sun
e_e = 0; % we assume circular orbit around the sun
a_e = 149597870.7; % [km] we use the value of 1au to set the semi-major axis of the earth's
% orbit

obliquity = deg2rad(23.44); % [rad] obliquity of the ecliptic
% shade of the earth is an ideal cylinder and only perturbation in 
% 2BP is solar radiation pressure

%% Part (a)

[r0,v0] = COE2rv(a,e,i,raan,omega,theta,muE);

[r0_ec,v0_ec] = EQ2EC(r0,v0,obliquity);

[a_ec,e_ec,i_ec,raan_ec,omega_ec,theta_ec] = rv2COE(muE,r0_ec,v0_ec);

i_ec = rad2deg(i_ec);
raan_ec = rad2deg(raan_ec);
omega_ec = rad2deg(omega_ec);
theta_ec = rad2deg(theta_ec);

%% Part (b)

% t = 0 when earth at summer solstice 

tau_e = 2*pi*sqrt(a_e^3/muS);
n = 2*pi / tau_e;
t = 0.3*tau_e;
F = srp(t,COE);

%%
% Plots
% Define the time range
t_values = linspace(0, tau_e, 10); % 1000 points from 0 to tau

% Initialize F_values for storing the results
F_values = zeros(3, length(t_values));

% Calculate F for each time value
for i = 1:length(t_values)
    F_values(:, i) = srp(t_values(i), COE); % CAMBIAR EL INPUT DE LOS COE!!!!!!!!!
end

% Plot the magnitude of F over time
figure;
plot(t_values, vecnorm(F_values));
title('Magnitude of F over time');
xlabel('Time [s]');
ylabel('Magnitude of F [N]');
figure;

% To illustrate the position of the earth
for i = 1:length(t_values)
    r_earth(:,i) = [-sin(2*pi*t_values(i)/tau_e), cos(2*pi*t_values(i)/tau_e), 0] * a_e;
    [r0, ~] = COE2rv(COE(1), COE(2), COE(3), COE(4), COE(5), COE(6) + n*t_values(i),muE);
    [r0_ec,~] = EQ2EC(r0,v0,obliquity);
    r_sat_sun(:,i) = r0_ec + r_earth(:,i);
    if norm(r_sat_sun(:,i)) > norm(r_earth(:,i))
        disp('SIII')
    end
end
plot(t_values,r_earth(1, :), 'r', t_values, r_earth(2, :), 'g', t_values, r_earth(3, :), 'b');
hold on
plot(t_values,r_sat_sun(1, :), t_values, r_sat_sun(2, :), t_values, r_sat_sun(3, :));

% Plot the components of F over time
% figure;
% plot(t_values, F_values(1, :), 'r');
% title('F_x over time');
% xlabel('Time [s]');
% ylabel('F [N]');
% legend('F_x')
% 
% figure;
% plot(t_values, F_values(2, :), 'r');
% title('F_y over time');
% xlabel('Time [s]');
% ylabel('F [N]');
% legend('F_y')
% 
% figure;
% plot(t_values, F_values(3, :), 'r');
% title('F_z over time');
% xlabel('Time [s]');
% ylabel('F [N]');
% legend('F_z')

figure;
plot(t_values, F_values(1, :), 'r', t_values, F_values(2, :), 'g', t_values, F_values(3, :), 'b');
title('Components of F over time');
xlabel('Time [s]');
ylabel('F [N]');
legend('$F_x$', '$F_y$', '$F_z$','Interpreter','Latex');