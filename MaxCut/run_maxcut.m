function R = run_maxcut(methods, nrep)
% RUN_MAXCUT  MaxCut problems of the rudy collection (130 instances of the
% BiqMac library): CD, ALM and PD from the solution of the SDP relaxation.
%
%   run_maxcut                     CD and ALM, followed by maxcut_table_paper,
%                                  plot_maxcut_paper and plot_maxcut_profiles
%   run_maxcut({'cd','alm','pd'})  any subset of the three methods
%   run_maxcut(methods, nrep)      running time = minimum over nrep runs
%                                  (default 3)
%
% PD needs about 130 s per instance (about 4.7 hours for the collection),
% so it is not run by default; the paper reports a single run of PD per
% instance, i.e. run_maxcut({'pd'}, 1).
%
% Per-instance results are saved after every instance in
% results/maxcut_rudy_<method>.mat; a method whose file is complete is not
% run again (delete the file to recompute), and an interrupted run is
% continued. Each method first makes one warm-up run on the first
% instance, which is not recorded. Every instance is solved nrep times;
% the runs are deterministic (this is checked), and the recorded running
% time is the minimum over the nrep runs. The starting points are computed
% once and cached (see maxcut_load).
if nargin < 1 || isempty(methods), methods = {'cd', 'alm'}; plots = true; else, plots = false; end
if nargin < 2 || isempty(nrep), nrep = 3; end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'common'));
resdir = fullfile(here, 'results');

% parameters of CD (all other parameters as in default_params)
par = default_params();
cd_tau = 10;  cd_eta = 0.5;  cd_c = 1;  cd_alpha = 0.3;  cd_beta0 = 5;

S = maxcut_load(par);
for m = 1:numel(methods)
	meth = lower(methods{m});
	assert(ismember(meth, {'cd', 'alm', 'pd'}), 'unknown method %s', meth);
	file = fullfile(resdir, ['maxcut_rudy_' meth '.mat']);
	if exist(file, 'file'), load(file, 'res'); else, res = struct([]); end
	p = par;
	if strcmp(meth, 'cd')
		p.tau = cd_tau;  p.eta = cd_eta;  p.alpha = cd_alpha;  p.beta0 = cd_beta0;
	end
	if numel(res) < numel(S)
		solve(meth, S(1), p, cd_c);      % warm-up run on the first instance
	end
	for i = numel(res)+1:numel(S)
		[x, o] = solve(meth, S(i), p, cd_c);
		for rep = 2:nrep
			[x2, o2] = solve(meth, S(i), p, cd_c);
			assert(isequal(x2, x) && o2.outer == o.outer && o2.inner == o.inner, ...
				'%s on %s is not deterministic', meth, S(i).name);
			o.time = min(o.time, o2.time);
		end
		r = rmfield(o, 'inner_flags');
		r.name = S(i).name;
		r.cut  = maxcut_cut(x, S(i).L);
		r.opt  = S(i).opt;
		if strcmp(meth, 'cd')
			r.par = struct('tau', cd_tau, 'eta', cd_eta, 'c', cd_c, 'alpha', cd_alpha, 'beta0', cd_beta0);
		end
		if isempty(res), res = r; else, res(i) = r; end
		save(file, 'res');
		fprintf('%-5s %-12s flag %d  outer %4d  inner %7d  cut %8g (%6.2f%%)  %8.2f s\n', ...
			meth, r.name, r.flag, r.outer, r.inner, r.cut, 100*r.cut/r.opt, r.time);
	end
	R.(meth) = res;
end

if plots
	maxcut_table_paper();
	plot_maxcut_paper();
	plot_maxcut_profiles();
end
end

function [x, o] = solve(meth, s, p, c)
% One run of method meth on instance s with parameters p.
switch meth
	case 'cd'
		[x, ~, o] = cdm(maxcut_problem(s.L, c), s.x0, p);
	case 'alm'
		[x, ~, o] = safeguarded_alm(maxcut_problem(s.L), s.x0, p);
	case 'pd'
		[x, ~, o] = pd(maxcut_problem(s.L), s.x0, p);
end
end
