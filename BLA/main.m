clc
clear
close all

training_ = false;

[training_data, test_data, Ts] = load_dataset();

nk = [1, 1];
if training_ == true
    orders = struc(3:7, 3:7);
    [nb, nf] = oestruc(training_data, test_data, orders, nk);
    model = oe(training_data, 'nb', nb, 'nf', nf, 'nk', nk);
else
    load model
    load nb
    load nf
end

compare(test_data, model, Inf)
