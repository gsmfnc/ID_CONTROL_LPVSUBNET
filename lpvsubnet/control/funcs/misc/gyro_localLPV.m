function [G] = gyro_localLPV(w1, q2, params, integrator)
    % GYRO_LOCALLPV Get linearized gyroscope model
    %
    %   Syntax:
    %       G = gyro_localLPV(w1, q2, params)
    %       G = gyro_localLPV(w1, q2, params, integrator)
    %
    %   Inputs:
    %       w1, q2: operating points
    %       params: parameters (from 'params.mat')
    %       integrator: output is angular position (true) or angular
    %           velocity (false)
    %
    %   Outputs:
    %       G: linear time-invariant model corresponding to selected
    %           operating points.
    %

    % Input validation
    narginchk(3, 4);
    if nargin <= 3
        integrator = true;
    end
    assert(islogical(integrator) && isscalar(integrator), ...
        '''integrator'' must be scalar logical');
    
    % q3 is fixed to 0 rad for the Sener D3 deliverable!
    q3 = 0;  % rad

    % Scheduling parameters
    p1 = cos(q2);
    p2 = sin(q2);
    p3 = cos(q3);
    p4 = sin(q3);
    p5 = w1;
    
    % Default Parameters
    I1 = params.Ic + params.Id;
    I2 = params.Jb + params.Jc + params.Jd;
    I3 = params.Id + params.Ka + params.Kb + params.Kc;
    I4 = params.Id - params.Jc - params.Jd + params.Kc;
    I5 = params.Ib + params.Ic - params.Kb - params.Kc;

    % Mass matrix
    M = [      params.Jd,      0, params.Jd*p1,                       params.Jd*p2*p3
                       0,     I1,            0,                                -I1*p4
            params.Jd*p1,      0, I4*p2^2 + I2,                          -I4*p1*p2*p3
         params.Jd*p2*p3, -I1*p4, -I4*p1*p2*p3, I4*p2^2*p4^2 - I4*p2^2 + I5*p4^2 + I3];
    
    % Coriolis matrix
    C = [params.fv1,                  0,                   0,                   0
                  0,         params.fv2,     params.Jd*p2*p5, -params.Jd*p1*p3*p5
                  0,   -params.Jd*p2*p5,          params.fv3,  params.Jd*p2*p4*p5
                  0, params.Jd*p1*p3*p5, -params.Jd*p2*p4*p5,          params.fv4];
    
    % State-space matrices
    Af = [zeros(4) eye(4); zeros(4) -M\C];
    Bf = [zeros(4); M\blkdiag(params.Km1, params.Km2, params.Km3, params.Km4)];
    
    % Remove states that are not required (e.g., red frame is fixed)
    A = Af([4 6 8], [4 6 8]);
    B = Bf([4 6 8], [2]);
    if integrator
        C = eye(1, 3);
    else
        C = [0, 0, 1];
    end
    D = zeros(1);
    
    % Linearized system
    G = ss(A, B, C, D);
    G.InputName = {'I2'};
    if integrator
        G.OutputName = {'q4'};
        G.StateName = {'q4', 'v2', 'v4'};
    else
        G.OutputName = {'v4'};
        G = minreal(G, [], false);
    end

end