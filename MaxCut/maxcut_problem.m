function prob = maxcut_problem(L, c)
% MAXCUT_PROBLEM  Rank-constrained reformulation of the MaxCut problem
% [Jia, Kanzow, Mehlitz, Wachsmuth, Math. Program. 2023, (6.3)],
%
%       min  -1/4*trace(L*W)
%       s.t. diag(W) = e,   W psd,   rank(W) <= 1,
%
% where L = Diag(A*e) - A is the Laplacian of the weighted graph. The
% variable W in S^N is stored as x = W(:) in R^(N^2). In the form
% min f(x) s.t. G(x) in C, x in D:
%
%       G(W) = diag(W) - e,   C = {0},   D = {W psd : rank(W) <= 1},
%
% with the projective mappings
%
%       Q_D(W)[Y] = c*(W*Y + Y*W),     Q_C(s) = 0.
%
% For W = lambda*v*v' in D with lambda > 0 and ||v|| = 1, the null space of
% Q_D(W) on S^N is {Y : Y*v = 0}, the normal space of D at W.

if nargin < 2 || isempty(c), c = 1; end

N  = size(L, 1);
L  = full((L + L')/2);
id = (0:N-1)'*N + (1:N)';                    % linear indices of diag(W) in W(:)
JG = sparse(1:N, id, 1, N, N^2);
gf = -L(:)/4;

prob.name = 'maxcut';
prob.n = N^2;
prob.m = N;
prob.N = N;
prob.c = c;
prob.L = L;

prob.f     = @(x) gf'*x;
prob.gradf = @(x) gf;

prob.G     = @(x) x(id) - 1;
prob.JG    = @(x) JG;
prob.hessG = [];                             % G is affine

prob.projD = @(x) proj_rank1_psd(x, N);
prob.projC = @(s) zeros(size(s));

P = sparse(reshape(reshape(1:N^2, N, N)', [], 1), 1:N^2, 1);   % P*vec(Y) = vec(Y')
prob.QD  = @(x, U) qd(x, U, N, c, P);
prob.dQD = @(x, u1, u2) dqd(u1, u2, N, c);
prob.QC  = @(s, W) zeros(N, size(W, 2));
prob.dQC = @(s, w1, w2) zeros(N, 1);
end

function V = qd(x, U, N, c, P)
% Q_D(W)*U, each column of U being vec(Y): vec(W*Y) for all columns at
% once, and vec(Y*W) = P*vec(W'*Y') with the commutation matrix P
W = reshape(x, N, N);
k = size(U, 2);
if issparse(U)
	% each nonzero u = Y(a,b) contributes u*W(:,a) to column b of W*Y
	% and u*W(b,:) to row a of Y*W (used for U = JG' with N nonzeros)
	[r, j, u] = find(U);
	a = mod(r-1, N) + 1;  b = (r - a)/N + 1;  q = numel(r);  t = (1:N)';
	I1 = t + (b' - 1)*N;              % entries of column b
	I2 = a' + (t - 1)*N;              % entries of row a
	V1 = W(:, a).*u';
	V2 = W(b, :)'.*u';
	J  = repmat(j', N, 1);
	V  = c*sparse([I1(:); I2(:)], [J(:); J(:)], [V1(:); V2(:)], N^2, k);
	return
end
V = reshape(W*reshape(U, N, N*k), N^2, k) ...
  + P*reshape(W'*reshape(P*U, N, N*k), N^2, k);
V = c*V;
end

function g = dqd(u1, u2, N, c)
% gradient w.r.t. W of <U1, Q_D(W)[U2]> = c*trace(U1'*(W*U2 + U2*W))
U1 = reshape(u1, N, N);  U2 = reshape(u2, N, N);
g  = c*reshape(U1*U2' + U2'*U1, [], 1);
end
