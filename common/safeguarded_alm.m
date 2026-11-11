function [x, lambda, out] = safeguarded_alm(prob, x0, par)
% SAFEGUARDED_ALM  Safeguarded augmented Lagrangian method of
% [Jia, Kanzow, Mehlitz, Wachsmuth, Math. Program. 2023, Algorithm 4.1] for
%
%       min f(x)   s.t.  G(x) in C,  x in D,
%
% with the augmented Lagrangian
%
%       L_rho(x,u) = f(x) + rho/2*dist_C(G(x) + u/rho)^2
%
% and the spectral gradient method as subproblem solver. The method stops
% when  V = ||G(x) - P_C(G(x) + u/rho)||_inf <= eps_out  [MP, (4.5)].
% The penalty parameter is multiplied by par.alm_tau whenever V fails to
% decrease by the factor par.alm_eta.
%
% Required fields of prob: n, m, f, gradf, G, JG, projD, projC.
%
% out.flag: 0 outer tolerance reached, 1 penalty bound exceeded,
%           2 evaluation budget exhausted, 3 max_outer reached

tstart = tic;

x  = prob.projD(x0);
Gx = prob.G(x);
rho = initial_penalty(prob.f(x), norm(Gx - prob.projC(Gx)));
rho0 = rho;
lambda = zeros(prob.m, 1);

k = 0;  inner = 0;  nfev = 1;  flag = 0;  V = Inf;  inner_flags = [];

if par.verbose
	fprintf('   k   inner  i_cum    f-ev.            f(x)       V_rho(x,u)         rho\n');
	fprintf('%4d %7s %6d %8d %15.8f %16s %11.3e\n', k, '---', inner, nfev, prob.f(x), '---', rho);
end

while true
	u   = proj_box(lambda, par.u_min, par.u_max);     % safeguarded multiplier
	tol = par.eps_inner/sqrt(k+1);
	phi = @(x) aug_lag(x, u, rho, prob);
	[x, sg] = spectral_gradient(phi, prob.projD, x, tol, par, par.max_fun - nfev);

	inner = inner + sg.iter;
	nfev  = nfev + sg.nfev;
	inner_flags(end+1) = sg.flag; %#ok<AGROW>
	k = k + 1;

	Gx     = prob.G(x);
	shift  = Gx + u/rho;
	Pshift = prob.projC(shift);
	lambda = rho*(shift - Pshift);
	V_old  = V;
	V      = norm(Gx - Pshift, inf);
	rho_used = rho;

	if par.verbose
		fprintf('%4d %7d %6d %8d %15.8f %16.3e %11.3e\n', k, sg.iter, inner, nfev, prob.f(x), V, rho);
	end

	if V <= par.eps_out,       flag = 0; break; end
	if nfev >= par.max_fun,    flag = 2; break; end
	if k >= par.max_outer,     flag = 3; break; end

	% penalty update (no update after the first outer iteration, V_old = Inf)
	if V > par.alm_eta*V_old
		rho = par.alm_tau*rho;
		if rho > par.beta_max, flag = 1; break; end
	end
end

out.flag   = flag;
out.outer  = k;
out.inner  = inner;
out.nfev   = nfev;
out.beta0  = rho0;
out.beta   = rho_used;                  % penalty parameter of the last subproblem
out.fval   = prob.f(x);
out.viol   = V;
out.infeas = norm(Gx - prob.projC(Gx), inf);   % dist_C(G(x)) in the inf-norm
out.inner_flags = inner_flags;
out.time   = toc(tstart);
end

function [L, gL] = aug_lag(x, u, rho, prob)
r = prob.G(x) + u/rho;
r = r - prob.projC(r);
L = prob.f(x) + 0.5*rho*(r'*r);
if nargout > 1
	gL = prob.gradf(x) + rho*(prob.JG(x)'*r);
end
end
