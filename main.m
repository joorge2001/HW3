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

% consider the orbit of the earth around the sun
e_e = 0; % we assume circular orbit around the sun
a_e = 149597870.7; % [km] we use the value of 1au to set the semi-major axis of the earth's
% orbit
tau_e = 2*pi*sqrt(a_e^3/muS);

obliquity = deg2rad(23.44); % [rad] obliquity of the ecliptic
% shade of the earth is an ideal cylinder and only perturbation in 
% 2BP is solar radiation pressure

%% Part (a)

[r0,v0] = COE2rv(a,e,i,raan,omega,theta,muE);

[r0_ec,v0_ec] = EQ2EC(r0,v0,obliquity);

[a_ec,e_ec,i_ec,raan_ec,omega_ec,theta_ec] = rv2COE(muE,r0_ec,v0_ec);

%% Part (b)

% F = srp(t);