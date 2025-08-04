function [Ap, Bp, Cp, Dp] = get_model(As, Bs, Cs, Ds)
    p1 = preal('p1', 'dt');
    p2 = preal('p2', 'dt');

    Ap = As{1} + As{2} * p1 + As{3} * p2;
    Bp = Bs{1} + Bs{2} * p1 + Bs{3} * p2;
    Cp = Cs{1} + Cs{2} * p1 + Cs{3} * p2;
    Dp = Ds{1} + Ds{2} * p1 + Ds{3} * p2;
end
