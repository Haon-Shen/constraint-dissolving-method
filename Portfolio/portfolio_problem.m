function prob = portfolio_problem(Sigma, mu, kappa, c)
% PORTFOLIO_PROBLEM  Sparse portfolio problem with an interval budget,
%
%       min  1/2*x'*Sigma*x - mu'*x
%       s.t. 0.8 <= e'*x <= 1,   x >= 0,   ||x||_0 <= kappa,
%
% in the form  min f(x)  s.t.  G(x) in C,  x in D  with
%
%       G(x) = e'*x - 1,   C = [-0.2, 0],   D = {x >= 0 : ||x||_0 <= kappa},
%
% and the projective mappings
%
%       Q_D(x) = c*diag(x_1^2,...,x_n^2),     Q_C(s) = (s + 0.2)*(-s).

if nargin < 4, c = 0.1; end

lo = -0.2;  up = 0;                  % C = [lo, up]
n  = numel(mu);
e  = ones(1, n);

prob.name  = 'portfolio';
prob.n = n;
prob.m = 1;
prob.c = c;
prob.kappa = kappa;
prob.budget = [1 + lo, 1 + up];      % bounds on e'*x

prob.f     = @(x) 0.5*(x'*(Sigma*x)) - mu'*x;
prob.gradf = @(x) Sigma*x - mu;

prob.G     = @(x) sum(x) - 1;
prob.JG    = @(x) e;
prob.hessG = [];                     % G is affine

prob.projD = @(x) proj_sparse_nonneg(x, kappa);
prob.projC = @(s) proj_box(s, lo, up);

prob.QD  = @(x, U) c*(x.^2).*U;
prob.dQD = @(x, u1, u2) 2*c*x.*u1.*u2;
prob.QC  = @(s, W) ((s - lo).*(up - s)).*W;
prob.dQC = @(s, w1, w2) (lo + up - 2*s).*w1.*w2;
end
