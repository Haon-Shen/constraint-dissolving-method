function [y, x, out] = pd(prob, x0, par)
% PD  Inexact penalty decomposition method of [Kanzow, Lapucci, Comput.
% Optim. Appl. 85 (2023), Algorithm 1] for
%
%       min f(x)   s.t.  G(x) in C,  x in D,
%
% based on the reformulation  min f(x) s.t. x - y = 0, G(x) in C, y in D,
% and the penalty function
%
%       q_beta(x,y) = f(x) + beta/2*( ||x-y||^2 + dist_C(G(x))^2 ).
%
% The subproblems are solved by the inexact alternating minimization of
% [Kanzow, Lapucci, Algorithm 2] with H = I: starting from (u,v) =
% (x^k,y^k), repeat
%
%       d = -grad_x q(u,v),   Armijo step u = u + t*d,   v = P_D(u)
%
% until ||grad_x q(u,v)||_inf <= delta_k, where delta_k =
% par.eps_inner/sqrt(k+1) is the inner tolerance of CD and ALM. The Armijo
% step is the largest t = par.pd_ls_beta^j, j = 0,1,..., with
% q(u+t*d,v) <= q(u,v) + par.pd_ls_gamma*t*g'*d. If the predicted decrease
% t*|g'*d| falls below the rounding level par.pd_ls_rtol*max(1,|q|), the
% subproblem is terminated.
%
% The penalty parameter is multiplied by par.pd_tau after every outer
% iteration. The method stops when
%
%       ||x - y||_inf + ||G(x) - P_C(G(x))||_inf <= eps_out.
%
% The starting point is x^0 = y^0 in P_D(x0), and the initial penalty
% parameter is given by initial_penalty.m (as for ALM). The computed
% solution is y (in D), returned first as by cdm and safeguarded_alm;
% out.infeas reports dist_C(G(y)), which is in general larger than the
% tolerance since the entries of x removed by P_D are lost.
%
% Required fields of prob: n, m, f, gradf, G, JG, projD, projC.
%
% out.flag: 0 outer tolerance reached, 1 penalty bound exceeded,
%           2 evaluation budget exhausted, 3 max_outer reached

tstart = tic;

y = prob.projD(x0);
x = y;
Gx   = prob.G(x);
beta = initial_penalty(prob.f(x), norm(Gx - prob.projC(Gx)));
beta0 = beta;

k = 0;  inner = 0;  nfev = 1;  flag = 0;  inner_flags = [];

if par.verbose
	fprintf('   k    alt  i_cum    f-ev.            f(y)     ||x-y||_inf   dist_C(G(x))        beta\n');
	fprintf('%4d %6s %6d %8d %15.8f %15s %14s %11.3e\n', k, '---', inner, nfev, prob.f(y), '---', '---', beta);
end

while true
	% subproblem: inexact alternating minimization (Algorithm 2)
	delta = par.eps_inner/sqrt(k+1);
	[x, y, lb] = alt_min(x, y, beta, delta, prob, par, par.max_fun - nfev);
	inner = inner + lb.iter;
	nfev  = nfev + lb.nfev;
	inner_flags(end+1) = lb.flag; %#ok<AGROW>
	k = k + 1;

	Gx = prob.G(x);
	dxy  = norm(x - y, inf);
	dC   = norm(Gx - prob.projC(Gx), inf);
	beta_used = beta;

	if par.verbose
		fprintf('%4d %6d %6d %8d %15.8f %15.3e %14.3e %11.3e\n', k, lb.iter, inner, nfev, prob.f(y), dxy, dC, beta);
	end

	if dxy + dC <= par.eps_out, flag = 0; break; end
	if nfev >= par.max_fun,     flag = 2; break; end
	if k >= par.max_outer,      flag = 3; break; end

	% penalty update
	beta = par.pd_tau*beta;
	if beta > par.beta_max, flag = 1; break; end
end

Gy = prob.G(y);
out.flag   = flag;
out.outer  = k;
out.alt    = inner;                     % number of alternations (= projections onto D)
out.inner  = inner;                     % number of x-steps (= alternations)
out.nfev   = nfev;
out.beta0  = beta0;
out.beta   = beta_used;                 % penalty parameter of the last subproblem
out.fval   = prob.f(y);
out.viol   = dxy + dC;                  % termination measure
out.dxy    = dxy;                       % ||x-y||_inf
out.infeas = norm(Gy - prob.projC(Gy), inf);   % dist_C(G(y)) in the inf-norm
out.inner_flags = inner_flags;          % flags of the subproblem solves (see alt_min)
out.time   = toc(tstart);
end


function [u, v, out] = alt_min(u, v, beta, delta, prob, par, max_fun)
% Algorithm 2 of [Kanzow, Lapucci] with H = I for min q_beta(x,y) s.t. y in D.
% out.flag: 0 ||grad_x q||_inf <= delta, 2 evaluation budget exhausted,
%           3 line search failed: the predicted decrease t*|g'*d| fell below
%             par.pd_ls_rtol*max(1,|q|), i.e., below the rounding level of q
[q, g] = pd_penalty(u, v, beta, prob);
nfev = 1;  iter = 0;  flag = 0;
while true
	if norm(g, inf) <= delta, flag = 0; break; end
	if nfev >= max_fun,       flag = 2; break; end
	d  = -g;
	gd = g'*d;
	t  = 1;
	qt = pd_penalty(u + t*d, v, beta, prob);  nfev = nfev + 1;
	ls_fail = false;
	while qt > q + par.pd_ls_gamma*t*gd
		t = par.pd_ls_beta*t;
		% the Armijo test cannot be decided in floating-point arithmetic
		if -t*gd <= par.pd_ls_rtol*max(1, abs(q)), ls_fail = true; break; end
		qt = pd_penalty(u + t*d, v, beta, prob);  nfev = nfev + 1;
	end
	if ls_fail, flag = 3; break; end
	u = u + t*d;
	v = prob.projD(u);
	[q, g] = pd_penalty(u, v, beta, prob);  nfev = nfev + 1;
	iter = iter + 1;
end
out.iter = iter;
out.nfev = nfev;
out.flag = flag;
end


function [q, gq] = pd_penalty(x, y, beta, prob)
% q_beta(x,y) and its gradient with respect to x.
Gx = prob.G(x);
r  = Gx - prob.projC(Gx);
dxy = x - y;
q = prob.f(x) + 0.5*beta*(dxy'*dxy + r'*r);
if nargout > 1
	gq = prob.gradf(x) + beta*(dxy + prob.JG(x)'*r);
end
end
