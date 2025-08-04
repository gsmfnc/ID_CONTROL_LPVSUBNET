function [data_train, data_test, Ts] = load_dataset()
    load("../data/ML_valid.mat");
    exp_p = [cos(exp_p(:, 2)), exp_p(:, 1), sin(exp_p(:, 2))];
    data_test = lpviddata(exp_y, exp_p, exp_u, Ts);

    load("../data/ML_estim.mat");
    exp_p = [cos(exp_p(:, 2)), exp_p(:, 1), sin(exp_p(:, 2))];
    data_train = lpviddata(exp_y, exp_p, exp_u, Ts);
end
