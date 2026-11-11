% RUN_ALL  Reproduce all tables and figures of the numerical section.
%
% Run this script from its own folder in MATLAB (R2025b; no toolbox is
% required). Each experiment writes its output to the subfolder results/
% of its folder; the outputs used in the paper are in the subfolders
% reference/. Approximate running times on a MacBook Air (Apple M4):
%
%   Rosenbrock MPCC (CD, ALM, PD)                      1 minute
%   portfolio problems with interval budget (CD, ALM)  0.5 minutes
%   mean-variance problems on market data (CD, ALM)    2 minutes (data on request)
%   MaxCut (CD, ALM; incl. the SDP starting points)    7 minutes
%   MaxCut (PD, see below)                             4.7 hours
%
% PD on the MaxCut problems is not run by default. To include it in the
% MaxCut table and bar chart, either set run_maxcut_pd = true below, or
% copy MaxCut/reference/maxcut_rudy_pd.mat to MaxCut/results/ before
% calling maxcut_table_paper and plot_maxcut_paper.

run_maxcut_pd = false;

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'common'));

% Interval-constrained Rosenbrock problem
cd(fullfile(root, 'MPCC'));
run_mpcc();
run_mpcc_pd();

% Sparse portfolio problems with an interval budget (orl and pard instances)
cd(fullfile(root, 'Portfolio'));
run_portfolio();

% Sparse mean-variance problems on market data
% (skipped if the stock return data are not available, see
% PortfolioMV/data/README.txt)
cd(fullfile(root, 'PortfolioMV'));
if isfile(fullfile('data', 'init_stock_data.mat'))
	run_mv();
else
	fprintf('PortfolioMV: stock return data not found, experiment skipped (see PortfolioMV/data/README.txt).\n');
end

% MaxCut problems (rudy collection)
cd(fullfile(root, 'MaxCut'));
run_maxcut({'cd', 'alm'});
if run_maxcut_pd
	run_maxcut({'pd'}, 1);
end
maxcut_table_paper();
plot_maxcut_paper();
plot_maxcut_profiles();

cd(root);
