function [best_order_nb, best_order_nf] = oestruc(training_data, test_data, ...
                                                  orders, nk)
    max_fit = 0;
    for i = 1:1:size(orders, 1)
        disp(i)
        for j = 1:1:size(orders, 1)
            nb = orders(i, :);
            nf = orders(j, :);

            model = oe(training_data, 'nb', nb, 'nf', nf, 'nk', nk);
            [yh, fit, x0] = compare(test_data, model, Inf);
            fits(i) = fit;
            if fit > max_fit
                max_fit = fit;
                best_order_nb = nb;
                best_order_nf = nf;
                disp(nb)
                disp(nf)
            end
        end
    end
end
