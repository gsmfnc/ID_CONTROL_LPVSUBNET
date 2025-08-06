close all

addpath('./funcs/lpvsubnet_scripts')

Ts = 0.01;
Ngrid = 20;

%% Non-normalized data
load('../../data/ML_estim.mat')
train_y_unnorm = exp_y;
train_u_unnorm = [exp_u exp_p(:, 1)];
load('../../data/ML_valid.mat')
test_y_unnorm = exp_y;
test_u_unnorm = [exp_u exp_p(:, 1)];

i2_norm_coeff = max(abs(train_u_unnorm(:, 1)));
q1dot_norm_coeff = max(abs(train_u_unnorm(:, 2)));
q4dot_norm_coeff = max(abs(train_y_unnorm));

clear exp_u exp_p exp_y

%% LPVcore results: bode diagram
% 
% load gyroscope_ss_model
% %load gyroscope_ss_model_noise
% sys_core = LPVcore.lpvss(ss_model.A, ss_model.B, ss_model.C, ss_model.D, Ts);
% 
% qd1_grid = linspace(30, 50, Ngrid);
% q2_grid = linspace(-pi/4, pi/4, Ngrid);
% 
% p1grid_core = cos(q2_grid);
% p2grid_core = qd1_grid;
% p3grid_core = sin(q2_grid);
% 
% fig1 = figure(1); sys_core.bode([p1grid_core; p2grid_core; p3grid_core]')
% %saveas(fig1, 'lpvcore_res', 'epsc')

%% LPVSUBNET results: bode diagram

if load_self_sched
    matrices
else
    matrices_ext
end

[Ap, Bp, Cp, Dp] = get_model(As, Bs, Cs, Ds);
encoder = {W0e, W1e, W2e, b0e, b1e, b2e, Ve, ce};
pnet = {W0, W1, W2, b0, b1, b2, V, c};
sys_m = {Ap, Bp, Cp, Dp};

[yks1, pout_max1, pout_min1] = ...
    compute_output(train_y, train_u, encoder, pnet, sys_m);
[yks2, pout_max2, pout_min2] = ...
    compute_output(valid_y, valid_u, encoder, pnet, sys_m);

%compare_output(train_y, yks1, 'fignum', 3)
%compare_output(valid_y, yks2, 'fignum', 4)

pout_max = max(pout_max1, pout_max2);
pout_min = min(pout_min1, pout_min2);

sys = LPVcore.lpvss(Ap, Bp, Cp, Dp, Ts);

% p1grid = linspace(pout_min(1), pout_max(1), Ngrid);
% p2grid = linspace(pout_min(2), pout_max(2), Ngrid);
% p3grid = linspace(pout_min(3), pout_max(3), Ngrid);

%fig2 = figure(2); sys.bode([p1grid; p2grid; p3grid]')
%saveas(fig2, 'lpvsubnet_res', 'epsc')

%% Changes by Chris from here:
% make new LPVSS system due to bug in LPVcore
p1 = preal('p1','dt','Range',[pout_min(1), pout_max(1)]);
p2 = preal('p2','dt','Range',[pout_min(2), pout_max(2)]);

Ac = Ap.matrices(:, :, 1) + Ap.matrices(:, :, 2) * p1 + ...
    Ap.matrices(:, :, 3) * p2;
Bc = Bp.matrices(:, :, 1) + Bp.matrices(:, :, 2) * p1 + ...
    Bp.matrices(:, :, 3) * p2;
Cc = Cp.matrices(:, :, 1) + Cp.matrices(:, :, 2) * p1 + ...
    Cp.matrices(:, :, 3) * p2;
Dc = Dp.matrices(:, :, 1) + Dp.matrices(:, :, 2) * p1 + ...
    Dp.matrices(:, :, 3) * p2;

G = LPVcore.lpvss(Ac, Bc, Cc, Dc, Ts);

% Add q4 as output and have a low-pass for q4dot to ensure parameter
% independent Cy
Gp = [c2d(zpk([], 0, 1), Ts); c2d(tf(100, [1 100]), Ts)] * G;

% Name inputs/outputs
Gp.InputName    = {'i2'};
Gp.OutputName   = {'q4', 'q4dot'};

% Scheduling-range
pRange = Gp.SchedulingTimeMap.Range;

global synthesize_flag
if synthesize_flag
    %% make generalized plant
    [P, ny, nu, MaxPole, AllWeights] = makeGeneralizedPlant(Gp, 1, 1);

    % Synthesis
    synopt = lpvsynOptions;
    synopt.solverOptions = sdpsettings('solver', 'mosek', 'verbose', 0);
    synopt.DirectFeedthrough = 0;

    disp('synopt.improveConditioning = 0;');
    synopt.improveConditioning = 0;

    [Kd, gam] = lpvsyn(P, ny, nu, synopt);

    %% Including scaling and filters in controller
    Ks = AllWeights.Mu*Kd*blkdiag(AllWeights.M, 1);

    % Rescale controller
    K = AllWeights.Su*Ks*(eye(ny)/AllWeights.Sy(1:ny, 1:ny));

    Ak = K.A;
    Bk = K.B;
    Ck = K.C;
    Dk = K.D;

    save ./funcs/mat_files/LPVcontroller K Ak Bk Ck Dk ny nu Ts pRange gam

else

    disp('synopt.improveConditioning = 0;');
    load ./funcs/mat_files/LPVcontroller

end
