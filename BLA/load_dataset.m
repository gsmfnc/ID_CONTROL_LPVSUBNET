function [training_data, test_data, Ts] = load_dataset()
    load("../data/ML_estim.mat")
    train_out = exp_y;
    train_in = [exp_u exp_p(:, 1)];

    load("../data/ML_valid.mat");
    test_out = exp_y;
    test_in = [exp_u exp_p(:, 1)];

    mm = mean(train_out);
    ss = std(train_out);
    train_out = (train_out - mm) / ss;
    test_out = (test_out - mm) / ss;

    mm = mean(train_in);
    ss = std(train_in);
    train_in(:, 1) = (train_in(:, 1) - mm(1)) / ss(1);
    train_in(:, 2) = (train_in(:, 2) - mm(2)) / ss(2);
    test_in(:, 1) = (test_in(:, 1) - mm(1)) ./ ss(1);
    test_in(:, 2) = (test_in(:, 2) - mm(2)) ./ ss(2);

    training_data = iddata(train_out, train_in, Ts);
    test_data = iddata(test_out, test_in, Ts);
end
