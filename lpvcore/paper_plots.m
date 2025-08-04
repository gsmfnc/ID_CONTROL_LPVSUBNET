close all
fsize = 18;
tsize = 12;
set(groot, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

fig = figure(1)
p1 = subplot(421);
plot(data_train.InputData); title('Training set', 'fontsize', tsize)
ylabel('$i_2$ (A)', 'interpreter', 'latex', 'fontsize', fsize)
set(gca, 'xticklabel', []);
yticks([-1, 0, 1, 2])
p2 = subplot(422);
plot(data_test.InputData); title('Test set', 'fontsize', tsize)
set(gca, 'xticklabel', []);
yticks([-1, 0, 1, 2])
p3 = subplot(423);
plot(data_train.SchedulingData(:, 1), 'g'); hold on
plot(data_train.SchedulingData(:, 3));
ylabel('[-]', 'interpreter', 'latex', 'fontsize', fsize)
legend('$\cos(q_2)$', '$\sin(q_2)$', 'interpreter', 'latex', 'fontsize', fsize)
set(gca, 'xticklabel', []);
yticks([-0.2, 0, 0.2, 0.4, 0.6, 0.8])
p4 = subplot(424);
plot(data_test.SchedulingData(:, 1), 'g'); hold on
plot(data_test.SchedulingData(:, 3));
legend('$\cos(q_2)$', '$\sin(q_2)$', 'interpreter', 'latex', 'fontsize', fsize)
set(gca, 'xticklabel', []);
yticks([-0.2, 0, 0.2, 0.4, 0.6, 0.8])
p5 = subplot(425);
plot(data_train.SchedulingData(:, 2));
ylabel('$\dot q_1$ (rad/s)', 'interpreter', 'latex', 'fontsize', fsize)
set(gca, 'xticklabel', []);
yticks([35, 40, 45])
p6 = subplot(426);
plot(data_test.SchedulingData(:, 2));
set(gca, 'xticklabel', []);
yticks([35, 40, 45])
p7 = subplot(427);
plot(data_train.OutputData);
ylabel('$\dot q_4$ (rad/s)', 'interpreter', 'latex', 'fontsize', fsize)
xlabel('$k$', 'interpreter', 'latex', 'fontsize', fsize)
p8 = subplot(428);
plot(data_test.OutputData);
xlabel('$k$', 'interpreter', 'latex', 'fontsize', fsize)
set(p1, 'position', [.1     .725 .4 .225])
set(p2, 'position', [.55    .725 .4 .225])
set(p3, 'position', [.1     .5 .4 .225])
set(p4, 'position', [.55    .5 .4 .225])
set(p5, 'position', [.1     .275 .4 .225])
set(p6, 'position', [.55    .275 .4 .225])
set(p7, 'position', [.1     .05 .4 .225])
set(p8, 'position', [.55    .05 .4 .225])

fig2 = figure(2)
[model_output, fit] = compare(data_test, oe_model, Inf);
plot(data_test.OutputData); hold on
plot(data_test.OutputData - model_output);
legend('$y$', '$y-\hat y$', 'interpreter', 'latex', 'fontsize', fsize)
bfr = compute_bfr(model_output, data_test.OutputData);
title("BFR=" + bfr + "%", 'interpreter', 'latex', 'fontsize', fsize)

fig3 = figure(3)
[model_output, fit] = compare(data_test, arx_model, Inf);
plot(data_test.OutputData); hold on
plot(data_test.OutputData - model_output);
legend('$y$', '$y-\hat y$', 'interpreter', 'latex', 'fontsize', fsize)
bfr = compute_bfr(model_output, data_test.OutputData);
title("BFR=" + bfr + "%", 'interpreter', 'latex', 'fontsize', fsize)

fig4 = figure(4)
[model_output, fit] = compare(data_test, ss_model, Inf);
plot(data_test.OutputData); hold on
plot(data_test.OutputData - model_output);
legend('$y$', '$y-\hat y$', 'interpreter', 'latex', 'fontsize', fsize)
bfr = compute_bfr(model_output, data_test.OutputData);
title("BFR=" + bfr + "%", 'interpreter', 'latex', 'fontsize', fsize)

function bfr = compute_bfr(y, yr)
num = 0;
den = 0;
mean_y = mean(yr);
num = norm(y - yr);
den = norm(yr - mean_y);

bfr = max(1 - num / den, 0) * 100
end
