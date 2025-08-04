load est_dat_1_noise.mat;

N=length(exp_y);
t = (0:Ts:(Ts*(N-1))).';
deg_scale=1/pi*180;
figure;
subplot(2,2,1);
plot(t,deg_scale*(out.q.Data(:,[2,4])));
title('Position')
legend('q2','q4');
xlabel('Time [s]');
ylabel('Amplitude [deg]');

subplot(2,2,2);
plot(t,(out.qd.Data(:,[2,4])));
title('Radial Speed')
legend('qd2','qd4');
xlabel('Time [s]');
ylabel('Speed [rad/s]');

subplot(2,2,3);
plot(t,[out.q1_ref.Data,out.qd.Data(:,1)]);
title('Radial Speed of The Flywheel')
legend('qd1 reference','qd1 actual');
xlabel('Time [s]');
ylabel('Speed [rad/s]');

subplot(2,2,4);
plot(t,exp_u);
title('Input for Gimbal 2')
xlabel('Time [s]');
ylabel('Current [A]');