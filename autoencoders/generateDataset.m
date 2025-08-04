clear
load("../data/ML_estim.mat")
Y = exp_y;
U = [exp_u exp_p(:, 1)];
load("../data/ML_valid.mat")
Y_val = exp_y;
U_val = [exp_u exp_p(:, 1)];
save('gyroscope.mat', 'U', 'Y', 'U_val', 'Y_val')
