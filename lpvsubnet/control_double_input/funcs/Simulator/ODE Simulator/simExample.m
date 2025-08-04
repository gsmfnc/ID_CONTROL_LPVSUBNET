clear
close all
clc

load('./../params.mat')

%% Initialization
Tend    = 10;               % End time
Tspan   = [0 Tend];         % Time span (can be a discrete time grid)
x0      = zeros(8, 1);      % Initial state

%% Solve ode45 script
% limitations:
% - Only supports constant inputs and disk references (no interpolation)
rDisk   = 40;               % Disk ref [rad/s]
u       = [0; 1; 0; 0];     % Actuate only the blue gimbal
lock    = [1 1 0 1];        % Lock the red gimbal; 1 = free, 0 = locked, [disk, bg, rg, sg]
[t, y]  = ode45(@(t, x) gyroNonlinSim(t, x, u, rDisk, params, lock), [0 Tend], x0);
