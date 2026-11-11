function prob = mv_problem(R, K, c)
% MV_PROBLEM  Sparse mean-variance portfolio problem with a return constraint,
%
%       min  1/2*x'*Sigma*x
%       s.t. mu'*x >= rho,   e'*x = 1,   x >= 0,   ||x||_0 <= K,
%
% where Sigma and mu are the sample covariance matrix and mean of the rows
% of the return matrix R, and rho = (sum(mu) - max(mu))/n. In the form
% min f(x) s.t. G(x) in C, x in D:
%
%       G(x) = (mu'*x - rho, e'*x - 1),   C = R_+ x {0},
%       D = {x >= 0 : ||x||_0 <= K},
%
% with the projective mappings
%
%       Q_D(x) = c*diag(x_1^2,...,x_n^2),     Q_C(s) = diag(s_1^2, 0).

if nargin < 3, c = 1; end

n     = size(R, 2);
mu    = mean(R, 1)';
Sigma = cov(R);
rho   = (sum(mu) - max(mu))/n;
JG    = [mu'; ones(1, n)];

prob.name = 'mv';
prob.n = n;
prob.m = 2;
prob.c = c;
prob.K = K;
prob.mu = mu;  prob.Sigma = Sigma;  prob.rho = rho;

prob.f     = @(x) 0.5*(x'*(Sigma*x));
prob.gradf = @(x) Sigma*x;

prob.G     = @(x) [mu'*x - rho; sum(x) - 1];
prob.JG    = @(x) JG;
prob.hessG = [];                     % G is affine

prob.projD = @(x) proj_sparse_nonneg(x, K);
prob.projC = @(s) [max(s(1), 0); 0];

prob.QD  = @(x, U) c*(x.^2).*U;
prob.dQD = @(x, u1, u2) 2*c*x.*u1.*u2;
prob.QC  = @(s, W) [s(1)^2*W(1,:); zeros(1, size(W, 2))];
prob.dQC = @(s, w1, w2) [2*s(1)*w1(1)*w2(1); 0];
end
