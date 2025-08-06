%% init
clear

addpath("./funcs")
load_gyro

load_self_sched = false; % load LPVSUBNET self-scheduling estimates
                         % NOTE: at the moment, externally-scheduled estimates
                         % do not work
synthesize_flag = false; % if false, it just loads previously synthesized
                         % controller
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
save('results_self.mat', 'extsched_gyro', 'selfsched_gyro');
% save('results_ext.mat', 'extsched_gyro', 'selfsched_gyro');

%% plotting
load results.mat
close all

plot_simres(extsched_gyro);
plot_simres(selfsched_gyro);

%% paper's plot
figure()
subplot(211)
plot(extsched_gyro.q4.Time, extsched_gyro.q4.Data); hold on
plot(extsched_gyro.q4ref.Time, extsched_gyro.q4ref.Data)
grid on
legend('$q_{4,ext}$ [rad/s]', '$q_{4,ref}$', 'interpreter', 'latex')

subplot(212)
plot(selfsched_gyro.q4.Time, selfsched_gyro.q4.Data); hold on
plot(selfsched_gyro.q4ref.Time, selfsched_gyro.q4ref.Data)
grid on
legend('$q_{4,self}$ [rad/s]', '$q_{4,ref}$', 'interpreter', 'latex')
xlabel('Time [s]')
