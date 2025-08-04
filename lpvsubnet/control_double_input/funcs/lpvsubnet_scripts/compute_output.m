function [yks, pout_max, pout_min] = compute_output(set_y, set_u, encoder, ...
                                                                    pnet, sys)
    Ap = sys{1}; Bp = sys{2}; Cp = sys{3}; Dp = sys{4};

    N = length(set_y);
    yks = zeros(N - 10, 1);
    pout_max = - Inf * ones(2, 1);
    pout_min = Inf * ones(2, 1);
    for i = 1:1:(N - 10)
        if i == 1
            iny = set_y(i:(i + 9));
            inu1 = set_u(i:(i + 9), 1);
            inu2 = set_u(i:(i + 9), 2);
            inu3 = set_u(i:(i + 9), 3);
            ein = [iny; inu1; inu2; inu3];
            eout = net_output(ein, encoder);
        else
            eout = xn;
        end

        x = eout;
        u = set_u(i + 9, :)';
        pin = [x; u];
        pout = [net_output(pin, pnet)]';
        pout_max = max(pout_max, pout');
        pout_min = min(pout_min, pout');

        yk = Cp.freeze(pout) * x + Dp.freeze(pout) * [u(1); u(2)];
        xn = Ap.freeze(pout) * x + Bp.freeze(pout) * [u(1); u(2)];

        yks(i) = yk;
    end
end
