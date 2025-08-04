function Kq2 = q2IdentController(fbw, params, int)
% Input parameters
%   fbw     Desired bandwidth [Hz]
%   params  Gyroscope parameter structure

    s = zpk('s');           % Laplace s

    % Model from i_2 -> q_2 (assuming other frames fixed)
    G = params.Km2/( (params.Ic + params.Id)*s^2 + params.fv2*s );
    
    %% Design low-frequent PD controller to keep q2 centered
    wbw     = 2*pi*fbw;     % Bandwidth frequency [rad/s]
    rwd     = 1/3;          % Bandwidth - lead ratio
    rwlp    = 6;            % Bandwidth - low-pass ratio
    blp     = 0.7;          % Low-pass normalized damping
    rwi     = 1/5;          % Bandwidth - integrator ratio    
    
    % Integrator
    if int
        ki      = 1/(sqrt(rwi^2 + 1));                                          % I-gain such that |Ki(wbw)| = 1
        Ki      = ki*(s + rwi*wbw)/s;
    else
        Ki      = 1;
    end
    
    % PD-lowpass controller
    kpd     = wbw*sqrt((2*blp*rwlp)^2+rwlp^4+1-2*rwlp^2)/sqrt(rwd^2+1);     % PD-gain such that |Kpd(wbw)| = 1
    Kpd     = kpd*(s + rwd*wbw)/(s^2 + 2*blp*rwlp*wbw*s + rwlp^2*wbw^2);

    % Kq2 controller
    kG      = 1/abs(evalfr(G, wbw));    % 1/|G(wbw)|
    Kq2     = kG*Ki*Kpd;
end