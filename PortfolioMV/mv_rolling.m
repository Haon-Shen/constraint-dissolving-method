function W = mv_rolling(ds, c, par, windows, nrep)
% MV_ROLLING  Rolling-window experiment for one data set ds (an entry of
% the output of mv_load_datasets).
%
% Window i (i = 1,...,T - T_est) uses the rows i,...,i+T_est-1 of ds.R as
% estimation data. The sparse mean-variance problem is solved by CD and
% ALM from the same starting point (solution of the problem without the
% cardinality constraint), and the resulting portfolios are evaluated on
% the return of the following period, row i+T_est.
%
%   windows - (optional) indices of the windows to be processed, default all
%   nrep    - (optional) repetitions for timing (minimum is reported), default 1
%
% W is a struct array over the processed windows.

[T, n] = size(ds.R);
N_out  = T - ds.T_est;
if nargin < 4 || isempty(windows), windows = 1:N_out; end
if nargin < 5, nrep = 1; end

W = struct('window', {}, 'rho', {}, 'x0', {}, 'x_cd', {}, 'x_alm', {}, ...
	'cd', {}, 'alm', {}, 'r_out', {});
xrel = [];
for j = 1:numel(windows)
	i    = windows(j);
	Rest = ds.R(i:i+ds.T_est-1, :);
	prob = mv_problem(Rest, ds.K, c);

	xrel = mv_warm_start(prob, par, xrel);       % warm-started by the previous window
	[xc, ~, oc] = cdm(prob, xrel, par);
	[xa, ~, oa] = safeguarded_alm(prob, xrel, par);
	for r = 2:nrep
		[~, ~, o] = cdm(prob, xrel, par);              oc.time = min(oc.time, o.time);
		[~, ~, o] = safeguarded_alm(prob, xrel, par);  oa.time = min(oa.time, o.time);
	end

	W(j).window = i;
	W(j).rho    = prob.rho;
	W(j).x0     = xrel;
	W(j).x_cd = xc;   W(j).x_alm = xa;
	W(j).cd     = add_info(oc, xc, prob);
	W(j).alm    = add_info(oa, xa, prob);
	W(j).r_out  = ds.R(i+ds.T_est, :);           % returns of the next period (percent)
end
end

function o = add_info(o, x, prob)
o.nnz        = nnz(x);
o.ret_slack  = prob.mu'*x - prob.rho;            % mu'*x - rho
o.budget     = sum(x);
o.inner_fail = any(o.inner_flags ~= 0);
end
