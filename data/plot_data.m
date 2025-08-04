estim_data = load('ML_estim.mat');
valid_data = load('ML_valid.mat');

Ts = 1/100;
t_estim = 0:Ts:((size(estim_data.exp_u, 1) - 1) * Ts);
t_valid = 0:Ts:((size(valid_data.exp_u, 1) - 1) * Ts);

t = tiledlayout(4, 2);

nexttile()
plot(t_estim, estim_data.exp_u)
ylabel('$i_2$ [A]', 'interpreter', 'latex')

nexttile()
plot(t_valid, valid_data.exp_u)

nexttile()
plot(t_estim, estim_data.exp_p(:, 1))
ylabel('$\dot q_1$ [rad/s]', 'interpreter', 'latex')

nexttile()
plot(t_valid, valid_data.exp_p(:, 1))

nexttile(); hold on
plot(t_estim, sin(estim_data.exp_p(:, 2)))
plot(t_estim, cos(estim_data.exp_p(:, 2)))
legend('$\sin(q_2)$', '$\cos(q_2)$', 'interpreter', 'latex')
ylim([0, 1.1])

nexttile(); hold on
plot(t_valid, sin(valid_data.exp_p(:, 2)))
plot(t_valid, cos(valid_data.exp_p(:, 2)))
legend('$\sin(q_2)$', '$\cos(q_2)$', 'interpreter', 'latex')
ylim([0, 1.1])

nexttile()
plot(t_estim, estim_data.exp_y)
xlabel('Time [s]')
ylabel('$\dot q_4$ [rad/s]', 'interpreter', 'latex')

nexttile()
plot(t_valid, valid_data.exp_y)
xlabel('Time [s]')

t.TileSpacing = 'compact';
