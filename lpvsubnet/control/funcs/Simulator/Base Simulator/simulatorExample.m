clear
close all
clc

%% Initialize model
load('params.mat')          % Load gyroscope parameters
load('Cdisk.mat')           % Load disk feedback controller

% Define initial conditions
% x0 = [q1 q2 q3 q4 w1 w2 w3 w4]
x0      = zeros(8, 1);      % Initialize states to 0
x0(5)   = 30;               % Initialize disk velocity to 30 [rad/s]

%% How to change parameters
% Look into the params struct and change a parameter, e.g., change
% disk friction:
params.fv1 = params.fv1*1.10;    % 10% increase

%% Simulate model
% Either open the simulation and manually run it or;
% Run from a script:
out = sim('GyroscopeSimulation.slx');

% state outputs
t = out.yout{1}.Values.Time;
x = out.yout{1}.Values.Data;