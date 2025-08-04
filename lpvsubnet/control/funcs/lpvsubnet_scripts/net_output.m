function out = net_output(in, net)
    W0 = net{1}; W1 = net{2}; W2 = net{3};
    b0 = net{4}; b1 = net{5}; b2 = net{6};
    V = net{7}; c = net{8};

    nonlin = W2 * tanh(W1 * tanh(W0 * in + b0) + b1) + b2;
    out = V * in + c + nonlin;
end
