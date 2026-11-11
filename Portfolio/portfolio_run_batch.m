function R = portfolio_run_batch(P, kappa, c, par, verbose)
% PORTFOLIO_RUN_BATCH  Run CD and ALM on all instances in P (output of
% portfolio_load_instances) with cardinality bound kappa and scaling c of
% Q_D. The running time of every run is the minimum over three repetitions
% (the runs are deterministic).
if nargin < 5, verbose = false; end
nrep = 3;
R = struct('name', {}, 'cd', {}, 'alm', {}, 'x_cd', {}, 'x_alm', {});

% warm-up runs, not timed
prob = portfolio_problem(P(1).Sigma, P(1).mu, kappa, c);
cdm(prob, P(1).x0, par);  safeguarded_alm(prob, P(1).x0, par);

for i = 1:numel(P)
	prob = portfolio_problem(P(i).Sigma, P(i).mu, kappa, c);
	[xc, ~, oc] = cdm(prob, P(i).x0, par);
	[xa, ~, oa] = safeguarded_alm(prob, P(i).x0, par);

	for r = 2:nrep
		[~, ~, o] = cdm(prob, P(i).x0, par);              oc.time = min(oc.time, o.time);
		[~, ~, o] = safeguarded_alm(prob, P(i).x0, par);  oa.time = min(oa.time, o.time);
	end

	R(i).name = P(i).name;
	R(i).cd = add_info(oc, xc);   R(i).alm = add_info(oa, xa);
	R(i).x_cd = xc;  R(i).x_alm = xa;
	if verbose
		fprintf('%-13s  CD f = %12.6f   ALM f = %12.6f\n', P(i).name, oc.fval, oa.fval);
	end
end
end

function o = add_info(o, x)
o.budget = sum(x);
o.nnz    = nnz(x);
o.inner_fail = any(o.inner_flags ~= 0);
end
