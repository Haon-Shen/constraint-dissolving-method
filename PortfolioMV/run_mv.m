function results = run_mv(names)
% RUN_MV  Sparse mean-variance portfolio selection with a return constraint
% on real stock data: CD versus the safeguarded ALM on all rolling windows,
% and out-of-sample evaluation of the computed portfolios.
%
%   run_mv           the five data sets of the paper: NASDAQ100, FTSE100 and
%                    the three S&P500 data sets (see mv_load_datasets)
%   run_mv(names)    selected data sets, e.g. run_mv({'nasdaq100'})
%
% Output files (folder results/):
%   mv_summary.txt    all tables in plain text
%   mv_tables.tex     rows of the LaTeX tables
%   fig_mv_out.{pdf,png}, fig_mv_regimes.{pdf,png}
%                     cumulative wealth curves (plot_mv_paper)
%   mv.mat            summary data (without the individual portfolios)

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'common'));
resdir = fullfile(here, 'results');
if ~exist(resdir, 'dir'), mkdir(resdir); end

if nargin < 1, names = {'nasdaq100', 'ftse100', 'sp_bear', 'sp_covid', 'sp_bull'}; end
DS = mv_load_datasets(names);

% problem-specific settings (all other parameters as in default_params)
par = default_params();
par.eps_out = 1e-5;
par.tau     = 10;                        % parameters of CD
par.eta     = 0.5;
par.alpha   = 0.01;
par.beta0   = 10;                        % (the formula of initial_penalty.m, used by
                                         % ALM, also gives 10 on all 1214 windows)
c           = 10;                        % scaling of Q_D

for d = 1:numel(DS)
	ds = DS(d);
	fprintf('%s: %d windows ...\n', ds.name, size(ds.R, 1) - ds.T_est);
	W = mv_rolling(ds, c, par);

	% timing of the first (in-sample) window: minimum over three repetitions
	W1 = mv_rolling(ds, c, par, 1, 3);
	W(1).cd.time = W1.cd.time;  W(1).alm.time = W1.alm.time;

	results(d).name = ds.name;  results(d).label = ds.label; %#ok<AGROW>
	results(d).n = size(ds.R, 2);  results(d).n_removed = ds.n_removed; %#ok<AGROW>
	results(d).T = size(ds.R, 1);  results(d).T_est = ds.T_est; %#ok<AGROW>
	results(d).K = ds.K;  results(d).NY = ds.NY;  results(d).c = c;  results(d).par = par; %#ok<AGROW>
	results(d).rho1 = W(1).rho;  results(d).nnz_relaxed1 = nnz(W(1).x0 > 1e-6); %#ok<AGROW>
	results(d).metrics = mv_metrics(ds, W); %#ok<AGROW>
end

save(fullfile(resdir, 'mv.mat'), 'results');
mv_write_output(results, resdir);
plot_mv_paper(resdir);
end
