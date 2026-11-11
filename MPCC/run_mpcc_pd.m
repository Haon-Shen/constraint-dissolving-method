function results = run_mpcc_pd()
% RUN_MPCC_PD  Penalty decomposition method (pd.m) on the interval-constrained
% Rosenbrock problem, from the starting points of run_mpcc. The results of
% CD and ALM are read from results/mpcc_rosenbrock.mat (output of run_mpcc,
% which must be run first) and are not recomputed. PD does not depend on
% the scaling c of Q_D.
%
% Output files (folder results/):
%   mpcc_rosenbrock_pd.txt   summary table with the rows of CD and ALM
%   mpcc_rosenbrock_pd.tex   row of PD in the format of mpcc_rosenbrock.tex
%   mpcc_rosenbrock_pd.mat   all per-run data of PD

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'common'));
resdir = fullfile(here, 'results');
B = load(fullfile(resdir, 'mpcc_rosenbrock.mat'));  B = B.results;

par = default_params();                  % eps_out = 1e-4 as in run_mpcc
fprintf('PD ...\n');
pd_res = mpcc_run_batch(@pd, mpcc_rosenbrock_problem(1), B.X0, par);

results.par = par;  results.seed = B.seed;  results.X0 = B.X0;  results.pd = pd_res;
save(fullfile(resdir, 'mpcc_rosenbrock_pd.mat'), 'results');

fid = fopen(fullfile(resdir, 'mpcc_rosenbrock_pd.txt'), 'w');
for out = [1, fid]
	fprintf(out, 'Interval-constrained Rosenbrock problem, %d starting points (seed %d)\n', size(B.X0, 2), B.seed);
	fprintf(out, 'PD: tau = %g, subproblems by [KL, Alg. 2] with H = I, gamma = %g, beta = %g, stop ||grad_x q||_inf <= %g/sqrt(k+1);  eps_out = %g on ||x-y||_inf + dist_C(G(x))\n', ...
		par.pd_tau, par.pd_ls_gamma, par.pd_ls_beta, par.eps_inner, par.eps_out);
	fprintf(out, 'CD and ALM from mpcc_rosenbrock.mat; for PD, Inner = alt. = x-steps (= projections onto D), infeas = dist_C(G(y))\n\n');
	fprintf(out, '%-6s %8s %9s %10s %10s %8s %7s %7s %8s %11s %11s %11s %9s\n', 'Alg.', 'c', 'Outer', ...
		'Inner', 'f-ev.', 'N_glob', 'N_loc', 'N_oth', 'N_succ', 'max infeas', 'avg infeas', 'time (s)', 'alt.');
	for j = 1:numel(B.c_list)
		print_row(out, 'CD', sprintf('%g', B.c_list(j)), B.cd(j));
	end
	print_row(out, 'ALM', '---', B.alm);
	print_row(out, 'PD', '---', pd_res);
	fprintf(out, '\nPD: runs with a subproblem not reaching the gradient tolerance: %d of %d\n', ...
		pd_res.n_inner_fail, numel(pd_res.flag));
	fprintf(out, 'PD: subproblems terminated by a failed line search (rounding level): %d of %d\n', ...
		pd_res.n_sub_lsfail, pd_res.n_sub);
end
fclose(fid);

fid = fopen(fullfile(resdir, 'mpcc_rosenbrock_pd.tex'), 'w');
fprintf(fid, '%% Alg. & c & Outer & Inner & N_global & N_local & mean time per run (ms) \\\\\n');
R = pd_res;
fprintf(fid, 'PD~\\cite{kanzow2023inexact} & -- & $%.3f$ & $%.3f$ & $%d$ & $%d$ & $%.2f$ \\\\\n', ...
	R.avg_outer, R.avg_inner, R.n_global, R.n_local, 1e3*R.avg_time);
fclose(fid);
end

function print_row(out, name, cstr, R)
fprintf(out, '%-6s %8s %9.3f %10.3f %10.1f %8d %7d %7d %8d %11.2e %11.2e %11.2e %9.1f\n', name, cstr, ...
	R.avg_outer, R.avg_inner, R.avg_nfev, R.n_global, R.n_local, R.n_other, ...
	R.n_success, R.max_infeas, R.avg_infeas, R.avg_time, R.avg_alt);
end
