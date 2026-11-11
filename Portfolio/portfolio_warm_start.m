function [x0, xrel, out] = portfolio_warm_start(prob, par)
% PORTFOLIO_WARM_START  Starting point for the sparse portfolio problem.
%
% The problem without the cardinality constraint,
%
%       min f(x)   s.t.  0.8 <= e'*x <= 1,  x >= 0,
%
% is a convex QP. It is solved to high accuracy by the spectral projected
% gradient method (no toolbox required); its solution xrel is returned
% unprojected in x0, since CD and ALM project the starting point onto D.

par.max_inner = 1e5;
tol  = 1e-8;
phi  = @(x) deal(prob.f(x), prob.gradf(x));
proj = @(x) proj_budget(x, prob.budget(1), prob.budget(2));

xinit = proj(ones(prob.n, 1)/prob.n);
[xrel, out] = spectral_gradient(phi, proj, xinit, tol, par);
x0 = xrel;
end
