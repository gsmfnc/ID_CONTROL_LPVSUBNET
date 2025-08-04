clc
clear
close all

value_for_ZeroIsNonFree = false;
plot_ = true;
training_ = false;

[data_train, data_test, Ts] = load_dataset();
plot_data(plot_, data_train, data_test)

na = 10;
nb = 10;
c2 = preal('cos(q2)', 'dt');
q1d = preal('q1d', 'dt');
s2 = preal('sin(q2)', 'dt');

%% LPV-ARX estimation
A = randn(1) + randn(1) * c2 + randn(1) * q1d + randn(1) * s2;
B = randn(1) + randn(1) * c2 + randn(1) * q1d + randn(1) * s2;
[A_poly, B_poly] = shift_pol(A, na, B, nb);

template_arx = lpvidpoly(A_poly, B_poly, [], [], [], 0, Ts, ...
    'ZeroIsNonFree', value_for_ZeroIsNonFree);
template_arx.InputName = {'Current Gimbal 2'};
template_arx.InputUnit = {'A'};
template_arx.OutputName = {'q4d'};
template_arx.OutputUnit = {'rad/s'};
if training_
    arx_model = lpvarx(data_train, template_arx);
else
    load gyroscope_arx_model
end
[rms_est, rms_test] = rms_computation(data_train, data_test, arx_model, plot_)

%% LPV-OE estimation
template_oe = lpvidpoly([], arx_model.B, [], [], arx_model.A, 0, Ts, ...
    'ZeroIsNonFree', value_for_ZeroIsNonFree);
template_oe.InputName = {'Current Gimbal 2'};
template_oe.InputUnit = {'A'};
template_oe.OutputName = {'q4d'};
template_oe.OutputUnit = {'rad/s'};
options_oe = lpvoeOptions('Display', 'off', ...
    'SearchOptions',  struct('StepSize', 10, 'StepTolerance', 1E-10), ...
    'Regularization', struct('Lambda', 1),...
    'Initialization', 'template');
if training_
    oe_model = lpvoe(data_train, template_oe, options_oe);
else
    load gyroscope_oe_model
end
[rms_est, rms_test] = rms_computation(data_train, data_test, oe_model, plot_)

%% LPV-PEM-SS estimation K=0
template_ss = lpvio2ss(arx_model.A, arx_model.B, na, nb, Ts, {c2, q1d, s2});
options_pem_ss = lpvssestOptions;
options_pem_ss.Display = 'on';
options_pem_ss.Initialization = 'template';
%options_pem_ss.SearchOptions.maxIter = 3;
if training_
    ss_model = lpvssest(data_train, template_ss, options_pem_ss);
else
    load gyroscope_ss_model
end
[rms_est, rms_test] = rms_computation(data_train, data_test, ss_model, plot_)

%%
save('gyroscope_arx_model', 'arx_model')
save('gyroscope_oe_model', 'oe_model')
save('gyroscope_ss_model', 'ss_model')
