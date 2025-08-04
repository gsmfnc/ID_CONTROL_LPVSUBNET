function [A_poly, B_poly] = shift_pol(A, na, B, nb);
    A_poly = eye(1);
    for k = 1:na
        if k == 1
            A_poly = {A_poly, pshift(A, -k)};
        else
            A_poly = {A_poly{:, 1:k}, pshift(A, -k)};
        end
    end
    B_poly = rand(1);
    for k = 1:nb
        if k == 1
            B_poly = {B_poly, pshift(B, -k)};
        else
            B_poly = {B_poly{:, 1:k}, pshift(B, -k)};
        end
    end
end
