function P = portfolio_load_instances(n, par, families)
% PORTFOLIO_LOAD_INSTANCES  Load the 20 instances of dimension n (families
% orl_05 and orl_005; letters a-j) together with their starting points,
% i.e. the solutions of the problems without the cardinality constraint.
% portfolio_load_instances(200, par, {'pard'}) loads the 10 pard instances.
if nargin < 3, families = {'orl_05', 'orl_005'}; end
letters  = 'abcdefghij';
P = struct('name', {}, 'Sigma', {}, 'mu', {}, 'x0', {});
for a = 1:numel(families)
	for b = 1:numel(letters)
		inst = read_portfolio_instance(n, families{a}, letters(b));
		prob = portfolio_problem(inst.Sigma, inst.mu, n);     % no cardinality constraint
		i = numel(P) + 1;
		P(i).name = inst.name;  P(i).Sigma = inst.Sigma;  P(i).mu = inst.mu;
		P(i).x0 = portfolio_warm_start(prob, par);
	end
end
end
