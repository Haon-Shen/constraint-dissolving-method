function S = mpcc_run_batch(solver, prob, X0, par)
% MPCC_RUN_BATCH  Run solver (@cdm, @safeguarded_alm or @pd) from every column
% of X0 and collect per-run data and summary statistics.

ntries  = size(X0, 2);
tol_sol = 1e-3;              % tolerance for identifying the limit point

S.X      = zeros(prob.n, ntries);
S.outer  = zeros(1, ntries);
S.inner  = zeros(1, ntries);
S.nfev   = zeros(1, ntries);
S.fval   = zeros(1, ntries);
S.infeas = zeros(1, ntries);
S.beta   = zeros(1, ntries);
S.time   = zeros(1, ntries);
S.flag   = zeros(1, ntries);
S.inner_fail = false(1, ntries);
S.alt    = zeros(1, ntries);  % alternations (pd only)
S.sub    = zeros(1, ntries);  % pd only: number of subproblems
S.sub_lsfail = zeros(1, ntries);  % pd only: subproblems ended by a failed line search
S.class  = zeros(1, ntries);  % 1 global minimizer, 2 local minimizer, 0 other

solver(prob, X0(:,1), par);            % warm-up run, not timed

for i = 1:ntries
	[x, ~, out] = solver(prob, X0(:,i), par);

	S.X(:,i)    = x;
	S.outer(i)  = out.outer;
	S.inner(i)  = out.inner;
	S.nfev(i)   = out.nfev;
	S.fval(i)   = out.fval;
	S.infeas(i) = out.infeas;
	S.beta(i)   = out.beta;
	S.time(i)   = out.time;
	S.flag(i)   = out.flag;
	S.inner_fail(i) = any(out.inner_flags ~= 0);
	if isfield(out, 'alt')
		S.alt(i) = out.alt;
		S.sub(i) = numel(out.inner_flags);
		S.sub_lsfail(i) = sum(out.inner_flags == 3);
	end

	if norm(x - prob.x_global, inf) <= tol_sol
		S.class(i) = 1;
	elseif norm(x - prob.x_local, inf) <= tol_sol
		S.class(i) = 2;
	end
end

S.avg_outer  = mean(S.outer);
S.avg_inner  = mean(S.inner);
S.avg_nfev   = mean(S.nfev);
S.avg_time   = mean(S.time);
S.n_global   = sum(S.class == 1);
S.n_local    = sum(S.class == 2);
S.n_other    = sum(S.class == 0);
S.n_success  = sum(S.flag == 0);     % runs satisfying the outer criterion
S.n_inner_fail = sum(S.inner_fail);  % runs with an unconverged inner solve
S.max_infeas = max(S.infeas);
S.avg_infeas = mean(S.infeas);
S.avg_alt    = mean(S.alt);
S.n_sub      = sum(S.sub);
S.n_sub_lsfail = sum(S.sub_lsfail);
end
