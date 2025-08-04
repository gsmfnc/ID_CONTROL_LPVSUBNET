function [params, Deviations] = getuncertainparams(params, worstcase, perturball)

if nargin < 3
    perturball = true;
end
if nargin < 2
    worstcase = false;
end

% sd = randi(50000); disp(sd)
sd = 35123;
rng(sd)

if perturball
    N = length(fieldnames(params));
else
    N = 7;
end

dev = (2*rand(N,1)-1);

if ~perturball
    % We have as uncertainties:
    %
    %       Parameter | I_B | I_C | I_D | fv1 | fv2 | fv4 | Km2 |
    % ----------------|-----|-----|-----|-----|-----|-----|-----|
    % Variation (+/-) | 20% | 20% | 20% | 20% | 20% | 20% | 10% |

    if worstcase
        deviations = [.2;.2;.2;.2;.2;.2;.1].*sign(dev);
    else
        deviations = [.4;.4;.4;.3;.6;.3;.3].*dev;
    end

    params.Ib = params.Ib*(1+deviations(1));
    params.Ic = params.Ic*(1+deviations(2));
    params.Id = params.Id*(1+deviations(3));
    params.fv1 = params.fv1*(1+deviations(4));
    params.fv2 = params.fv2*(1+deviations(5));
    params.fv4 = params.fv4*(1+deviations(6));
    params.Km2 = params.Km2*(1+deviations(7));

else
    % we assume ALL parameters are +- 40% deviated
    paramnames = fieldnames(params);
    deviations = zeros(N,1);
    for ii = 1:N
        deviations(ii) = 0.5*sign(dev(ii));
        params.(paramnames{ii}) = params.(paramnames{ii})*(1+deviations(ii));
    end
end

if nargout > 1
    Deviations = deviations;
end

end