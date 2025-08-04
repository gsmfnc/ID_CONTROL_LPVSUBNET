%% "single" input
tmp = load('ML_valid.mat');
load('ML_estim_noise_free.mat');

tmp.Noise.snr

exp_p_1_noise = randn(size(exp_p(:, 1))) * 3e-2;
snr(exp_p(:, 1), exp_p_1_noise)
exp_p(:, 1) = exp_p(:, 1) + exp_p_1_noise;

exp_p_2_noise = randn(size(exp_p(:, 2))) * 2e-3;
snr(exp_p(:, 2), exp_p_2_noise)
exp_p(:, 2) = exp_p(:, 2) + exp_p_2_noise;

exp_y_noise = randn(size(exp_y)) * 5e-3;
snr(exp_y, exp_y_noise)
exp_y = exp_y + exp_y_noise;

save('ML_estim.mat', 'exp_p', 'exp_y', 'exp_u')

%% "double" input
load('double_input_data.mat');

exp_p_1_noise = randn(size(exp_p(:, 1))) * 0.75;
snr(exp_p(:, 1), exp_p_1_noise)
exp_p(:, 1) = exp_p(:, 1) + exp_p_1_noise;

exp_p_2_noise = randn(size(exp_p(:, 2))) * 2e-3;
snr(exp_p(:, 2), exp_p_2_noise)
exp_p(:, 2) = exp_p(:, 2) + exp_p_2_noise;

exp_y_noise = randn(size(exp_y)) * 5e-3;
snr(exp_y, exp_y_noise)
exp_y = exp_y + exp_y_noise;

exp_p_est = exp_p(1:30000, :);
exp_y_est = exp_y(1:30000, :);
exp_u_est = exp_u(1:30000, :);
exp_p_val = exp_p(30001:60000, :);
exp_y_val = exp_y(30001:60000, :);
exp_u_val = exp_u(30001:60000, :);

save('DBL_estim.mat', 'exp_p_est', 'exp_y_est', 'exp_u_est')
save('DBL_valid.mat', 'exp_p_val', 'exp_y_val', 'exp_u_val')
