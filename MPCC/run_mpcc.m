function results = run_mpcc(ntries, c_list)
% RUN_MPCC  Interval-constrained Rosenbrock problem (an MPCC): CD for
% several scalings c of Q_D, and the safeguarded ALM, from the same
% starting points. The results of PD are computed by run_mpcc_pd.
%
%   run_mpcc                  1000 starting points, c in {0.01,0.1,1,10,100}
%   run_mpcc(ntries, c_list)  custom setting
%
% Output files (folder results/):
%   mpcc_rosenbrock.txt   summary table
%   mpcc_rosenbrock.tex   rows of the LaTeX table (CD for each c, and ALM)
%   mpcc_rosenbrock.mat   all per-run data

if nargin < 1 || isempty(ntries), ntries = 1000; end
if nargin < 2 || isempty(c_list), c_list = [0.01, 0.1, 1, 10, 100]; end

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'common'));
resdir = fullfile(here, 'results');
if ~exist(resdir, 'dir'), mkdir(resdir); end

% parameters of CD (all other parameters as in default_params)
par  = default_params();
par.tau   = 10;
par.eta   = 0.5;
par.alpha = 10;
seed = 2021;
X0   = mpcc_starting_points(ntries, seed);
base = 'mpcc_rosenbrock';

% ----- CD for all values of c -----
for j = 1:numel(c_list)
	prob = mpcc_rosenbrock_problem(c_list(j));
	fprintf('CD, c = %g ...\n', c_list(j));
	cd_res(j) = mpcc_run_batch(@cdm, prob, X0, par); %#ok<AGROW>
end

% ----- ALM with its default parameters (does not depend on c) -----
fprintf('ALM ...\n');
alm_res = mpcc_run_batch(@safeguarded_alm, mpcc_rosenbrock_problem(1), X0, par);

results.par = par;  results.seed = seed;  results.X0 = X0;
results.c_list = c_list;  results.cd = cd_res;  results.alm = alm_res;
save(fullfile(resdir, [base '.mat']), 'results');

% ----- summary -----
fid = fopen(fullfile(resdir, [base '.txt']), 'w');
for out = [1, fid]
	fprintf(out, 'Interval-constrained Rosenbrock problem, %d starting points (seed %d)\n', ntries, seed);
	fprintf(out, 'CD: tau = %g, eta = %g, alpha = %g;  ALM: tau = %g, eta = %g;  eps_out = %g, eps_inner = %g/sqrt(k+1)\n\n', ...
		par.tau, par.eta, par.alpha, par.alm_tau, par.alm_eta, par.eps_out, par.eps_inner);
	fprintf(out, '%-6s %8s %9s %10s %10s %8s %7s %7s %8s %11s %11s\n', 'Alg.', 'c', 'Outer', ...
		'Inner', 'f-ev.', 'N_glob', 'N_loc', 'N_oth', 'N_succ', 'max infeas', 'time (s)');
	for j = 1:numel(c_list)
		print_row(out, 'CD', sprintf('%g', c_list(j)), cd_res(j));
	end
	print_row(out, 'ALM', '---', alm_res);
end
fclose(fid);

fid = fopen(fullfile(resdir, [base '.tex']), 'w');
fprintf(fid, '%% Alg. & c & Outer & Inner & N_global & N_local & mean time per run (ms) \\\\\n');
for j = 1:numel(c_list)
	tex_row(fid, 'CD', sprintf('$%g$', c_list(j)), cd_res(j));
end
fprintf(fid, '\\midrule\n');
tex_row(fid, 'ALM~\cite{jia2023augmented}', '--', alm_res);
fclose(fid);
end

function print_row(out, name, cstr, R)
fprintf(out, '%-6s %8s %9.3f %10.3f %10.1f %8d %7d %7d %8d %11.2e %11.2e\n', name, cstr, ...
	R.avg_outer, R.avg_inner, R.avg_nfev, R.n_global, R.n_local, R.n_other, ...
	R.n_success, R.max_infeas, R.avg_time);
end

function tex_row(fid, name, cstr, R)
fprintf(fid, '%s & %s & $%.3f$ & $%.3f$ & $%d$ & $%d$ & $%.2f$ \\\\\n', name, cstr, ...
	R.avg_outer, R.avg_inner, R.n_global, R.n_local, 1e3*R.avg_time);
end
