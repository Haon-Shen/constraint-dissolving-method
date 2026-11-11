function [x0, out] = mv_warm_start(prob, par, xinit)
% MV_WARM_START  Starting point for the sparse mean-variance problem: the
% solution of the problem without the cardinality constraint,
%
%       min 1/2*x'*Sigma*x   s.t.  mu'*x >= rho,  e'*x = 1,  x >= 0,
%
% which is a convex QP. It is solved to high accuracy by the safeguarded ALM
% with D = R^n_+ (no toolbox required). CD and ALM project x0 onto D.
%
% xinit (optional) is the starting point of this auxiliary solve.

if nargin < 3 || isempty(xinit), xinit = ones(prob.n, 1)/prob.n; end

relaxed = prob;
relaxed.projD = @(x) max(x, 0);

par.eps_out   = 1e-8;
par.eps_inner = 1e-6;
par.alm_tau   = 10;
par.alm_eta   = 0.8;
par.verbose   = 0;
[x0, ~, out] = safeguarded_alm(relaxed, xinit, par);
end
