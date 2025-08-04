function [P, ny, nu, MaxPole, AllWeights] = makeGeneralizedPlant(Gp, isDiscrete, addLowPass)
% Function to construct generalized plant. This function also contains all
% the design choices. 

%% Preallocs
if nargin < 2
    isDiscrete = false;
end
if nargin < 3
    addLowPass = false;
end
AllWeights = {};

%% Settings
% all the settings for the generalized plant are configured here
umax    = 10;
ymax    = [45*(pi/180), 3]; % q4 dev and q4dot dev
wB      = 1*(2*pi);

if isDiscrete
    c2d_options = c2dOptions(...
        'Method', 'zoh', ... %options: zoh, foh, impulse, tustin, matched, least-squares
        'PrewarpFrequency', 0, ...
        'FitOrder','auto', ... used for least-squares method
        'FractDelayApproxOrder', 0);
end

if addLowPass
    wLP = 15*wB;
    Mu = tf(wLP,[1 wLP]);
    if isDiscrete
        Mu = c2d(Mu,Gp.Ts);   % Don't use tustin as results in non-zero D
    end
else
    Mu = 1;
end

%% Scaling
Su      = diag(umax);
Sy      = diag(ymax);
G       = (Sy\eye(size(Sy)))*Gp*Su; 

%% Generalized plant
WeFull = db2mag(-5)*tf([1 wB],[1 0]);   % WeFull := We*M
[~,Mi,We] = rncf(WeFull);
M = inv(Mi,'min');

Wq4dot = 1;

Wr = 1;
Wd = 0.1;

Wu   = makeweight(db2mag(0),[10*wB,3],db2mag(20));

if isDiscrete
    WeFull = c2d(WeFull, Gp.Ts, c2d_options);
    M = c2d(M, Gp.Ts, c2d_options);
    We = c2d(We, Gp.Ts, c2d_options);
    Wu = c2d(Wu, Gp.Ts, c2d_options);
end

ny = 2;
nu = 1;

M.InputName = {'e'};
M.OutputName = {'y'};

G.InputName = Gp.InputName;
G.OutputName = Gp.OutputName;

% systemnames = 'G M';                        %#ok
% inputvar = '[r(1);d(1);u(1)]';              %#ok
% outputvar = '[M; G(2:2); u; M; G(2:2)]';    %#ok
% input_to_G = '[u+d]';                       %#ok
% input_to_M = '[r-G(1:1)]';                  %#ok
% cleanupsysic = 'yes';                       %#ok
Puw = connect(G, M, ...
    sumblk('i2 = u + d'), ...
    sumblk('e = r - q4'), ...
    {'r', 'd', 'u'}, ...
    {'y', 'q4dot', 'u', 'y', 'q4dot'});


Wz = blkdiag(We, Wq4dot, Wu);
Ww = blkdiag(Wr, Wd);

Wzy = blkdiag(Wz, eye(ny));
Wwu = blkdiag(Ww, eye(nu) * Mu);

P = Wzy * Puw * Wwu;

AllWeights.Su = Su;
AllWeights.Sy = Sy;
AllWeights.WeFull = WeFull;
AllWeights.M = M;
AllWeights.Mu = Mu;
AllWeights.We = We;
AllWeights.Wq4dot = Wq4dot;
AllWeights.Wr = Wr;
AllWeights.Wd = Wd;
AllWeights.Wu = Wu;
AllWeights.Wz = Wz;
AllWeights.Ww = Ww;
AllWeights.Wzy = Wzy;
AllWeights.Wwu = Wwu;

if isDiscrete
    MaxPole = NaN;
else
    MaxPole = abs(eig(Wu))+1;
end

end
