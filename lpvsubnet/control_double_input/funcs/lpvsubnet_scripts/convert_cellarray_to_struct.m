function netstruct = convert_cellarray_to_struct(net)
    netstruct.W0 = net{1}; netstruct.W1 = net{2}; netstruct.W2 = net{3};
    netstruct.b0 = net{4}; netstruct.b1 = net{5}; netstruct.b2 = net{6};
    netstruct.V = net{7};  netstruct.c = net{8};
end
