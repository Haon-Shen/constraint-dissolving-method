function [x, s, out] = cdm(prob, x0, par)
% CDM  Constraint dissolving method CD of the paper for
%
%       min f(x)   s.t.  G(x) in C,  x in D.
%
% The subproblems  min h_beta(x,s)  s.t. (x,s) in D x C  are solved by the
% spectral gradient method. The starting point is x^0 in P_D(x0),
% s^0 = P_C(G(x^0)). The penalty parameter is multiplied by par.tau
% whenever ||G(x)-s|| fails to decrease by the factor par.eta.
%
% Besides the fields listed in cdm_penalty, prob provides the projections
% projD(x) and projC(s). The parameters par are described in default_params.
%
% out.flag: 0 outer tolerance reached, 1 penalty bound exceeded,
%           2 evaluation budget exhausted, 3 max_outer reached

tstart = tic;
n = prob.n;

x = prob.projD(x0);
s = prob.projC(prob.G(x));
projK = @(z) [prob.projD(z(1:n)); prob.projC(z(n+1:end))];

H     = prob.G(x) - s;
viol  = norm(H);
if isfield(par, 'beta0') && ~isempty(par.beta0)
	beta = par.beta0;                   % fixed initial penalty parameter
else
	beta = initial_penalty(prob.f(x), viol);
end
beta0 = beta;

k = 0;  inner = 0;  nfev = 1;  flag = 0;  inner_flags = [];

if par.verbose
	fprintf('   k   inner  i_cum    f-ev.            f(x)   ||G(x)-s||_inf        beta\n');
	fprintf('%4d %7s %6d %8d %15.8f %16.3e %11.3e\n', k, '---', inner, nfev, prob.f(x), norm(H,inf), beta);
end

while true
	% approximate M-stationary point of the penalty subproblem
	tol = par.eps_inner/sqrt(k+1);
	phi = @(z) cdm_penalty(z, beta, par.alpha, prob);
	[z, sg] = spectral_gradient(phi, projK, [x; s], tol, par, par.max_fun - nfev);
	x = z(1:n);  s = z(n+1:end);

	inner = inner + sg.iter;
	nfev  = nfev + sg.nfev;
	inner_flags(end+1) = sg.flag; %#ok<AGROW>
	k = k + 1;

	H        = prob.G(x) - s;
	viol_old = viol;
	viol     = norm(H);
	beta_used = beta;

	if par.verbose
		fprintf('%4d %7d %6d %8d %15.8f %16.3e %11.3e\n', k, sg.iter, inner, nfev, prob.f(x), norm(H,inf), beta);
	end

	if norm(H, inf) <= par.eps_out,  flag = 0; break; end
	if nfev >= par.max_fun,          flag = 2; break; end
	if k >= par.max_outer,           flag = 3; break; end

	% penalty update
	if viol > par.eta*viol_old
		beta = par.tau*beta;
		if beta > par.beta_max,      flag = 1; break; end
	end
end

out.flag   = flag;
out.outer  = k;
out.inner  = inner;
out.nfev   = nfev;
out.beta0  = beta0;
out.beta   = beta_used;                 % penalty parameter of the last subproblem
out.fval   = prob.f(x);
out.viol   = norm(H, inf);              % ||G(x)-s||_inf
Gx = prob.G(x);
out.infeas = norm(Gx - prob.projC(Gx), inf);   % dist_C(G(x)) in the inf-norm
out.inner_flags = inner_flags;
out.time   = toc(tstart);
end
