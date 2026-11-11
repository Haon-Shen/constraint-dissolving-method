function results = run_portfolio(n_list, kappa_list, family)
% RUN_PORTFOLIO  Sparse portfolio problems with an interval budget: CD
% versus the safeguarded ALM on the Frangioni-Gentile instances
% (portfolio_problem.m).
%
%   run_portfolio                        everything for the paper: the orl and
%                                        pard instances for n in {200,300,400}
%                                        and kappa in {3,5,7}, followed by
%                                        portfolio_table_paper (table) and
%                                        plot_portfolio_paper (figure)
%   run_portfolio(n_list, kappa_list, family)
%                                        one family, 'orl' (20 instances per n:
%                                        orl<n>_05_a-j and orl<n>_005_a-j) or
%                                        'pard' (10 instances per n: pard<n>_a-j)
%
% Both methods start from the projection onto D of the solution of the
% problem without the cardinality constraint (portfolio_warm_start). ALM
% uses its default parameters (beta0 by the formula of initial_penalty.m,
% which gives 10 on all these problems). The running time of every run is
% the minimum over three repetitions (portfolio_run_batch).
%
% Output files (folder results/<family>/):
%   portfolio_n<n>_kappa<kappa>.txt   per-instance results and summary
%   portfolio.mat                     all data

if nargin == 0                           % everything for the paper
	results = run_portfolio([], [], 'orl');
	run_portfolio([], [], 'pard');
	portfolio_table_paper();
	plot_portfolio_paper();
	return
end
if nargin < 1 || isempty(n_list), n_list = [200, 300, 400]; end
if nargin < 2 || isempty(kappa_list), kappa_list = [3, 5, 7]; end
if nargin < 3 || isempty(family), family = 'orl'; end

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'common'));
resdir = fullfile(here, 'results', family);
if ~exist(resdir, 'dir'), mkdir(resdir); end

% problem-specific settings (all other parameters as in default_params)
par = default_params();
par.eps_out = 1e-5;
par.tau   = 2;                           % parameters of CD
par.eta   = 0.5;
par.alpha = 0.001;
par.beta0 = 300;
c         = 0.01;                        % scaling of Q_D

switch family
	case 'orl',  families = {'orl_05', 'orl_005'};
	case 'pard', families = {'pard'};
	otherwise,   error('unknown family %s', family);
end

j = 0;
for n = n_list
	P = portfolio_load_instances(n, par, families);
	for kappa = kappa_list
		R = portfolio_run_batch(P, kappa, c, par, true);

		j = j + 1;
		results(j).n = n;  results(j).kappa = kappa;  results(j).c = c; %#ok<AGROW>
		results(j).par = par;  results(j).runs = R; %#ok<AGROW>
		results(j).summary = portfolio_summary(R); %#ok<AGROW>

		fid = fopen(fullfile(resdir, sprintf('portfolio_n%d_kappa%d.txt', n, kappa)), 'w');
		for out = [1, fid]
			print_table(out, n, kappa, c, par, R, results(j).summary);
		end
		fclose(fid);
	end
end

save(fullfile(resdir, 'portfolio.mat'), 'results');
end

function print_table(out, n, kappa, c, par, R, s)
fprintf(out, '\nSparse portfolio with interval budget: n = %d, kappa = %d, c = %g\n', n, kappa, c);
fprintf(out, 'CD: tau = %g, eta = %g, alpha = %g, beta0 = %g;  ALM: tau = %g, eta = %g;  eps_out = %g, eps_inner = %g/sqrt(k+1)\n\n', ...
	par.tau, par.eta, par.alpha, par.beta0, par.alm_tau, par.alm_eta, par.eps_out, par.eps_inner);
fprintf(out, '%-13s | %12s %8s %4s %6s %7s %9s %4s | %12s %8s %4s %6s %7s %9s %4s\n', 'instance', ...
	'f CD', 'e''x', 'nnz', 'outer', 'inner', 'time (s)', 'flag', ...
	'f ALM', 'e''x', 'nnz', 'outer', 'inner', 'time (s)', 'flag');
for i = 1:numel(R)
	C = R(i).cd;  A = R(i).alm;
	fprintf(out, '%-13s | %12.6f %8.5f %4d %6d %7d %9.4f %4d | %12.6f %8.5f %4d %6d %7d %9.4f %4d\n', R(i).name, ...
		C.fval, C.budget, C.nnz, C.outer, C.inner, C.time, C.flag, ...
		A.fval, A.budget, A.nnz, A.outer, A.inner, A.time, A.flag);
end
fprintf(out, '\nlower objective value: CD %d, ALM %d, tie %d (relative difference <= %g)\n', ...
	s.cd_wins, s.alm_wins, s.ties, s.tol_tie);
fprintf(out, 'mean objective value: CD %.4f, ALM %.4f\n', s.mean_f);
fprintf(out, 'mean outer / inner / f-ev. / time: CD %.2f / %.1f / %.1f / %.4f s,  ALM %.2f / %.1f / %.1f / %.4f s\n', ...
	s.outer(1), s.inner(1), s.nfev(1), s.time(1), s.outer(2), s.inner(2), s.nfev(2), s.time(2));
fprintf(out, 'runs meeting the outer criterion: CD %d, ALM %d of %d;  runs with unconverged inner solves: CD %d, ALM %d\n', ...
	s.success, numel(R), s.inner_fail);
end
