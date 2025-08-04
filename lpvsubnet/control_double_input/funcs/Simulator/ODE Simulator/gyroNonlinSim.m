function [xdot] = gyroNonlinSim(t, x, u, rDisk, pars, lock)
    % Assign states
    q2 = x(2);      % Position blue gimbal
    q3 = x(3);      % Position red gimbal
    w1 = x(5);      % Velocity disk
    w2 = x(6);      % Velocity blue gimbal
    w3 = x(7);      % Velocity red gimbal
    w4 = x(8);      % Velocity silver gimbal
    qdot = [w1; w2; w3; w4];
    u = u(:);
    
    % Extract what frames are locked
    if nargin == 5
        idxL = find(lock);
    else
        idxL = 1:4;
    end

    % Extract parameters from structure
    cellfun(@(x) assignin('caller', x, pars.(x)), fieldnames(pars));
    
    % Mass matrix
    M(1, 1) = Jd;
    M(1, 2) = 0;
    M(1, 3) = Jd*cos(q2);
    M(1, 4) = Jd*cos(q3)*sin(q2);
    M(2, 1) = M(1, 2);
    M(2, 2) = Ic + Id;
    M(2, 3) = 0;
    M(2, 4) = -sin(q3)*(Ic + Id);
    M(3, 1) = M(1, 3);
    M(3, 2) = M(2, 3);
    M(3, 3) = Jb + Jc + Jd + Id*sin(q2)^2 - Jc*sin(q2)^2 - Jd*sin(q2)^2 + Kc*sin(q2)^2;
    M(3, 4) = -cos(q2)*cos(q3)*sin(q2)*(Id - Jc - Jd + Kc);
    M(4, 1) = M(1, 4);
    M(4, 2) = M(2, 4);
    M(4, 3) = M(3, 4);
    M(4, 4) = Id + Ka + Kb + Kc + Ib*sin(q3)^2 + Ic*sin(q3)^2 - Id*sin(q2)^2 + Jc*sin(q2)^2 + Jd*sin(q2)^2 - Kb*sin(q3)^2 - Kc*sin(q2)^2 - Kc*sin(q3)^2 + Id*sin(q2)^2*sin(q3)^2 - Jc*sin(q2)^2*sin(q3)^2 - Jd*sin(q2)^2*sin(q3)^2 + Kc*sin(q2)^2*sin(q3)^2;
    
    % Coriolis matrix
    C(1, 1) = 0;
    C(1, 2) = (Jd*w4*cos(q2)*cos(q3))/2 - (Jd*w3*sin(q2))/2;
    C(1, 3) = -(Jd*sin(q2)*(w2 + w4*sin(q3)))/2;
    C(1, 4) = (Jd*w2*cos(q2)*cos(q3))/2 - (Jd*w3*sin(q2)*sin(q3))/2;
    C(2, 1) = (Jd*w3*sin(q2))/2 - (Jd*w4*cos(q2)*cos(q3))/2;
    C(2, 2) = 0;
    C(2, 3) = (Jd*w1*sin(q2))/2 - w4*((cos(q3)*(Ic + Id))/2 - (cos(q2)^2*cos(q3)*(Id - Jc - Jd + Kc))/2 + (cos(q3)*sin(q2)^2*(Id - Jc - Jd + Kc))/2) - w3*cos(q2)*sin(q2)*(Id - Jc - Jd + Kc);
    C(2, 4) = w4*cos(q2)*cos(q3)^2*sin(q2)*(Id - Jc - Jd + Kc) - (Jd*w1*cos(q2)*cos(q3))/2 - w3*((cos(q3)*(Ic + Id))/2 - (cos(q2)^2*cos(q3)*(Id - Jc - Jd + Kc))/2 + (cos(q3)*sin(q2)^2*(Id - Jc - Jd + Kc))/2);
    C(3, 1) = -(Jd*sin(q2)*(w2 - w4*sin(q3)))/2;
    C(3, 2) = w4*((cos(q3)*(Ic + Id))/2 - (cos(q2)^2*cos(q3)*(Id - Jc - Jd + Kc))/2 + (cos(q3)*sin(q2)^2*(Id - Jc - Jd + Kc))/2) - (Jd*w1*sin(q2))/2 + w3*cos(q2)*sin(q2)*(Id - Jc - Jd + Kc);
    C(3, 3) = (w2*sin(2*q2)*(Id - Jc - Jd + Kc))/2;
    C(3, 4) = w2*((cos(q3)*(Ic + Id))/2 - (cos(q2)^2*cos(q3)*(Id - Jc - Jd + Kc))/2 + (cos(q3)*sin(q2)^2*(Id - Jc - Jd + Kc))/2) - w4*cos(q3)*sin(q3)*(Ib + Ic - Kb - Kc + Id*sin(q2)^2 - Jc*sin(q2)^2 - Jd*sin(q2)^2 + Kc*sin(q2)^2) + (Jd*w1*sin(q2)*sin(q3))/2;
    C(4, 1) = (Jd*w2*cos(q2)*cos(q3))/2 - (Jd*w3*sin(q2)*sin(q3))/2;
    C(4, 2) = (Jd*w1*cos(q2)*cos(q3))/2 - w3*((cos(q3)*(Ic + Id))/2 + (cos(q2)^2*cos(q3)*(Id - Jc - Jd + Kc))/2 - (cos(q3)*sin(q2)^2*(Id - Jc - Jd + Kc))/2) - w4*cos(q2)*cos(q3)^2*sin(q2)*(Id - Jc - Jd + Kc);
    C(4, 3) = w4*cos(q3)*sin(q3)*(Ib + Ic - Kb - Kc + Id*sin(q2)^2 - Jc*sin(q2)^2 - Jd*sin(q2)^2 + Kc*sin(q2)^2) - w2*((cos(q3)*(Ic + Id))/2 + (cos(q2)^2*cos(q3)*(Id - Jc - Jd + Kc))/2 - (cos(q3)*sin(q2)^2*(Id - Jc - Jd + Kc))/2) - (Jd*w1*sin(q2)*sin(q3))/2 + w3*cos(q2)*sin(q2)*sin(q3)*(Id - Jc - Jd + Kc);
    C(4, 4) = w3*cos(q3)*sin(q3)*(Ib + Ic - Kb - Kc + Id*sin(q2)^2 - Jc*sin(q2)^2 - Jd*sin(q2)^2 + Kc*sin(q2)^2) - w2*cos(q2)*cos(q3)^2*sin(q2)*(Id - Jc - Jd + Kc);
    
    % Viscous frictions
    D = blkdiag(fv1, fv2, fv3, fv4);

    % Input matrix (motor constants)
    Km = blkdiag(Km1, Km2, Km3, Km4);

    % Eliminate states that are locked
    Ml = M(idxL, idxL);
    Cl = C(idxL, idxL);
    Dl = D(idxL, idxL);
    Kml = Km(idxL, idxL);
    
    % Control disk such that it tracks the reference "rDisk"
    pDisk = 0.3;
    u(1) = -pDisk*(w1-rDisk) + rDisk*(fv1/Km1);

    % State equations
    xdot = zeros(8, 1);
    xdot(1:4, 1) = qdot;
    xdot(idxL+4, 1) = -Ml\(Cl+Dl)*qdot(idxL) + Ml\(Kml*u(idxL));
end

