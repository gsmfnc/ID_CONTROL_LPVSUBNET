%% init
clear

addpath("./funcs")
load_gyro

synthesize_flag = true;
synthesis_main

x0 = [0; 0; 0; 0; 50; 0; 0; 0];

type = 0; % external-scheduling
type = 1; % self-scheduling

encoder_struct = convert_cellarray_to_struct(encoder);
pnet_struct = convert_cellarray_to_struct(pnet);

clc

tfin = 50;

%% simulations
gyro_sim = 'gyroscope.slx';
load_system(gyro_sim);

type = 0; % external-scheduling
extsched_gyro = sim(gyro_sim, tfin);

type = 1; % self-scheduling
selfsched_gyro = sim(gyro_sim, tfin);

%% save results
save('results.mat', 'extsched_gyro', 'selfsched_gyro');

%% plotting
load results.mat
close all

plot_simres(extsched_gyro);
plot_simres(selfsched_gyro);
