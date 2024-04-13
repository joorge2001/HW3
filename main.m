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
S = 20; % [m^2] always points towards the sun
% orbit around earth

rE = 6371; % [km]
muE = 3.986e5; % [km3/s2]
muS = 1.32712440018e11; % [km3/s2]

COE = hw3data();

name = {'Earth';'Satellite'};
title_graph = {'Initial situation in an ECI-eq RF'};
ellipticOrbitPlotter(COE,muE,name,title_graph,rE,1)

% consider the orbit of the earth around the sun
e_e = 0; % we assume circular orbit around the sun
a_e = 149597870.7; % [km] we use the value of 1au to set the semi-major axis of the earth's
% orbit

obliquity = deg2rad(23.44); % [rad] obliquity of the ecliptic
% shade of the earth is an ideal cylinder and only perturbation in 
% 2BP is solar radiation pressure

%% Part (a)

[r0,v0] = COE2rv(COE(1),COE(2),COE(3),COE(4),COE(5),COE(6),muE);

[r1_ec,v1_ec] = EQ2EC(r0,v0,obliquity);

[a_ec,e_ec,i_ec,raan_ec,omega_ec,theta_ec] = rv2COE(muE,r1_ec,v1_ec);

COE2 = [a_ec,e_ec,i_ec,raan_ec,omega_ec,theta_ec];
name = {'Earth';'Satellite ec'};
title_graph = {'Initial situation in an ECI-ec RF'};
ellipticOrbitPlotter(COE2,muE,name,title_graph,rE,1)


i_ec = rad2deg(i_ec);
raan_ec = rad2deg(raan_ec);
omega_ec = rad2deg(omega_ec);
theta_ec = rad2deg(theta_ec);

%% Part (b)

% t = 0 when earth at summer solstice 
tau = 2*pi*sqrt(COE(1)^3/muE);
tau_e = 2*pi*sqrt(a_e^3/muS);
n = 2*pi / tau;
n_e = 2*pi / tau_e;

%%
% Plots
% Define the time range
t_values = linspace(0, tau_e, 1000); % 1000 points from 0 to tau

% Initialize F_values for storing the results
F_values = zeros(3, length(t_values));

% Calculate F for each time value
for i = 1:length(t_values)
    F_values(:, i) = srp(t_values(i));
end

% Plot the magnitude of F over time
figure;
plot(t_values, vecnorm(F_values));
title('Magnitude of F over time');
xlabel('Time [s]');
ylabel('Magnitude of F [N]');

% To illustrate the position of the earth and satellite
for i = 1:length(t_values)
    r_earth(:,i) = [-sin(2*pi*t_values(i)/tau_e), cos(2*pi*t_values(i)/tau_e), 0] * a_e;
    [r0, v0] = COE2rv(COE(1), COE(2), COE(3), COE(4), COE(5), COE(6) + n*t_values(i),muE);
    r0_plot(:,i) = r0;
    [r1_ec(:,i),~] = EQ2EC(r0,v0,obliquity);
    r_sat_sun(:,i) = r1_ec(:,i) + r_earth(:,i);
end

figure;
plot3(r_earth(1, :), r_earth(2, :),r_earth(3, :), 'b');
hold on
plot3(r_sat_sun(1, :), r_sat_sun(2, :),r_sat_sun(3, :), 'r');
title('Position of the Earth and the satellite that orbits it over time ');
legend('Earth','Satellite');
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]')

figure;
plot(t_values, r_sat_sun(1, :), 'r', t_values, r_sat_sun(2, :), 'g', t_values, r_sat_sun(3, :), 'b');
title('Components of $\overrightarrow{r}^{s}_{2}$ over time');
xlabel('Time [s]');
ylabel('r [km]');
legend('$r_x$', '$r_y$', '$r_z$');

figure;
plot(t_values, F_values(1, :), 'r', t_values, F_values(2, :), 'g', t_values, F_values(3, :), 'b');
title('Components of F over time');
xlabel('Time [s]');
ylabel('F [N]');
legend('$F_x$', '$F_y$', '$F_z$','Interpreter','Latex');