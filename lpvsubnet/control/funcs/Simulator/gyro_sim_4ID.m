%% Initialisation
clearvars; close all; clc;
addpath(genpath('ID Simulator'));
GyroscopeSimulationFile = 'GyroscopeSimulationID.slx';
load('params.mat');
load('Cdisk.mat');
Ts = 1E-2;  % Sampling Time - Fized step solver is used
flag=[];
Noise=[];
x0=zeros(8,1);
%% Experiment Desing
directory = 'data/';
namedataset = 'val';
savenow = false;
switch namedataset
    case 'val'
        rng(12); 
        nsets = 20; 
    case 'train'
        rng(13); 
        nsets = 200; 
end

N = 5000;  % Number of data samples (simulation steps) / period;
NumPeriod = 1;

u_mag=1;                % input magnitude
q1_mag=40;              % q1_ref nominal 

flag.noise = true;      % Simulate with output nosie
Noise.Variance = 7*1E-7; % Output nosie variance (in avg. 35dB SNR)

flag.input = 'SRGS';    %'RGS': Generates a Random, Gaussian Signal.
                        %'RBS': Generates Random, Binary Signal.
                        %'SINE': Generates a sum-of-sinusoid signal.
                        %'SRGS': Sine carrier with random Gaussian signal

flag.q1_ref= 'STCASE';  %'STCASE': staircase signal with multiple leves
                        %'CONST': constant signal level 

x0(5)=40;               % start system with fylwheel spinning   
x0(2)=0/180*pi;         % to start with an intial q2 position          

SNRs = [];
for datasetid=1:nsets

%% Data Generation

switch flag.input
    case 'RGS'
        u_sigma = u_mag;
        u_Band = [0,1];
        u_Range= [-u_sigma,u_sigma];
        u=idinput([N, 1, NumPeriod],flag.input,u_Band,u_Range);
    case 'RBS'
        u_Band = [0,3/50];
        u_Range= [-u_mag,u_mag];
        u=idinput([N, 1, NumPeriod],flag.input,u_Band,u_Range);
    case 'SINE'
        u_Band = [1/N,30/50];
        u_Range= [-u_mag,u_mag];
        NumSinusoids = 20;
        NumTrials = 20;
        GridSkip = 1;
        SineData = [NumSinusoids,NumTrials,GridSkip];
        u=idinput([N, 1, NumPeriod],'sine',u_Band,u_Range,SineData);
        disp(['Frequency band: [' num2str(u_Band/(2*Ts)) ']Hz']);
    case 'SRGS'
        freq_base=(1+rand(1))*2*pi; % 1-2Hz carrier signal
        u_base=(u_mag/2)*sin(freq_base*0:Ts:(Ts*(N-1)))';
        u_sigma = u_mag/3;
        u_Band = [0,1];
        u_Range= [-u_sigma,u_sigma];
        u=u_base+idinput([N, 1, NumPeriod],'RGS',u_Band,u_Range);
end
N=N*NumPeriod;
t = (0:Ts:(Ts*(N-1))).';
u = timeseries(u,t);

switch flag.q1_ref
    case 'CONST'
        q1_ref= timeseries(q1_mag*ones(N,1),t);
    case 'STCASE'
        steady_min =200*2;
        steady_max =400*2;
        q1_Range=[30,50];
        k=0;
        steady_length=steady_max-steady_min;
        range_length=q1_Range(2)-q1_Range(1);
        q1_ref=[];
        while k<N
            k_delta=steady_min+round(rand(1)*steady_length);
            q1_ref=[q1_ref;(rand(1)*range_length+q1_Range(1))*ones(k_delta,1)];
            k=k+k_delta;
        end
        q1_ref= timeseries(q1_ref(1:N,1),t);
end
options = simset('SrcWorkspace', 'current');
out = sim(['ID Simulator/', GyroscopeSimulationFile], [], options);
    
q = out.q.Data;
qd = out.qd.Data;

exp_u=u.Data;
exp_p=[qd(:,1),q(:,2)];
if flag.noise
    %Noise.Variance=((norm(qd(:,4),2))/db2mag(35)/sqrt(N))^2;
    Noise.e=randn(N,1)*sqrt(Noise.Variance);
    exp_y=qd(:,4)+Noise.e;
    Noise.snr=snr(qd(:,4),Noise.e);
    SNRs = [SNRs; Noise.snr];
else
    exp_y=qd(:,4);
end
%IDdat=iddata(exp_y,exp_u,Ts);

%save('est_dat_1_noise.mat','exp_y','exp_u','exp_p','Ts','out');

if savenow
    filename = [directory namedataset '-' num2str(datasetid) '.mat'];
    %save('../../../data-nocontrol-short/val_dat_9.mat','exp_y','exp_u','exp_p','Ts','out', 'q1_ref');
    save(filename,'exp_y', 'exp_u', 'exp_p','Ts', 'out');
end

end
figure; plot(SNRs)

%% Plotting
deg_scale=1/pi*180;
figure;
subplot(2,2,1);
plot(t,deg_scale*(q(:,[2,4])));
title('Position')
legend('q2','q4');
xlabel('Time [s]');
ylabel('Amplitude [deg]');

subplot(2,2,2);
plot(t,(qd(:,[2,4])));
title('Radial Speed')
legend('qd2','qd4');
xlabel('Time [s]');
ylabel('Speed [rad/s]');

subplot(2,2,3);
plot(t,[q1_ref.Data,qd(:,1)]);
title('Radial Speed of The Flywheel')
legend('qd1 reference','qd1 actual');
xlabel('Time [s]');
ylabel('Speed [rad/s]');

subplot(2,2,4);
plot(t,u.Data);
title('Input for Gimbal 2')
xlabel('Time [s]');
ylabel('Current [A]');
