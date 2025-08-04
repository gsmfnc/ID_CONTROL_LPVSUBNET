function plot_data(plot_, data_train, data_test)
    if plot_
        figure(1)
        subplot(211); plot(data_train.OutputData)
        title('Output training data')
        subplot(212); plot(data_test.OutputData)
        title('Output test data')

        figure(2)
        subplot(321); plot(data_train.SchedulingData(:, 1))
        title('Scheduling 1 training data')
        subplot(322); plot(data_test.SchedulingData(:, 1))
        title('Scheduling 1 test data')
        subplot(323); plot(data_train.SchedulingData(:, 2))
        title('Scheduling 2 training data')
        subplot(324); plot(data_test.SchedulingData(:, 2))
        title('Scheduling 2 test data')
        subplot(325); plot(data_train.SchedulingData(:, 3))
        title('Scheduling 3 training data')
        subplot(326); plot(data_test.SchedulingData(:, 3))
        title('Scheduling 3 test data')

        figure(3)
        subplot(211); plot(data_train.InputData)
        title('Input training data')
        subplot(212); plot(data_test.InputData)
        title('Input test data')
    end
end
