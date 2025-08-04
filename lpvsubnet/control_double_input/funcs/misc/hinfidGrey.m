function [sys_id, gamma, info] = hinfidGrey(G, p, Minit, Delta, opts)
% HINFID Computes lpvlfr model given N LTI systems and scheduling values
% p.
%
%   Syntax:
%
%       [sys_id, gamma, info] = hinfid(G, p, sys_init, opts);
%       [sys_id, gamma, info] = hinfid(G, p, sys_init);
%
%   Inputs:
%       G: cell array of LTI models (tf, ss, ...).
%       p: cell array of frozen scheduling variables.
%       sys_init: initial LPV-LFR estimate.
%       opts: options for HINFSTRUCT. Default: same as HINFSTRUCT defaults.
%
%   Outputs:
%       sys_id: estimated LPV-LFR model.
%       gamma: maximum error between given LTI systems G and sys_id.
%       info: information structure returned by HINFSTRUCT.

% Validate
% if nargin == 3
%     opts = hinfstructOptions;
% end
% if length(G) ~= length(p)
%     error('Number of input models must equal number of frozen scheduling variables');
% end
% for i=1:numel(p)
%     if numel(p{i}) ~= sys_init.Np
%         error(['Each supplied frozen scheduling variable must have ', int2str(sys_init.Np), ...
%             ' elements.']);
%     end
% end

% Get dimensions
nx = size(Minit.A, 1);
ny = 1;
nu = 1;
nD = size(Minit.D, 1) - ny;

% create matrix with unknown parameters and loop over all models 
Ghat = [];
for i=1:length(G)   
    deltaP = freeze(Delta,p{i});    
    Gp = G{i}-lft(tf(deltaP),Minit); 
    Ghat = blkdiag(Ghat,Gp);
end

% solve using hinfstruct
[Gstar, gamma, info] = hinfstruct(Ghat, opts);
A = Gstar.Blocks.A.Value;
Bw = Gstar.Blocks.Bw.Value;
Bu = Gstar.Blocks.Bu.Value;
Cz = Gstar.Blocks.Cz.Value;
Cy = Gstar.Blocks.Cy.Value;
Dzw = zeros(nD);
Dzu = zeros(nD, nu);
Dyw = zeros(ny, nD);
Dyu = zeros(ny, nu);

sys_id = lpvlfr(Delta, ss(A, [Bw Bu], [Cz; Cy], [Dzw Dzu; Dyw Dyu]));

end