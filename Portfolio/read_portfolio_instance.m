function inst = read_portfolio_instance(n, family, letter)
% READ_PORTFOLIO_INSTANCE  Read one instance of the Frangioni-Gentile
% portfolio collection from the folder data/ (see data/README.txt).
%
%   n      - 200, 300 or 400
%   family - 'orl_05', 'orl_005' or 'pard'
%   letter - 'a', ..., 'j'
%
% Returns inst.name, inst.Sigma (covariance matrix) and inst.mu (expected
% returns). The files .rho and .bds of the collection (required return
% and bounds on the entries) are not used by the interval-budget model.

datadir = fullfile(fileparts(mfilename('fullpath')), 'data');

switch family
	case 'orl_05',  name = sprintf('orl%d_05_%s',  n, letter);
	case 'orl_005', name = sprintf('orl%d_005_%s', n, letter);
	case 'pard',    name = sprintf('pard%d_%s',    n, letter);
	otherwise,      error('unknown family');
end
if n == 200 && ~strcmp(family, 'pard')
	name = strrep(name, '_', '-');      % file naming of the collection
end

% .mat (text file): dimension followed by the entries of Sigma
d = scanfile(fullfile(datadir, [name '.mat']));
assert(d(1) == n && numel(d) == n^2 + 1, 'unexpected format of %s.mat', name);
Sigma = reshape(d(2:end), n, n);
Sigma = (Sigma + Sigma')/2;

% .txt: dimension followed by n pairs, the first entry of each is mu_i
d = scanfile(fullfile(datadir, [name '.txt']));
assert(d(1) == n && numel(d) == 2*n + 1, 'unexpected format of %s.txt', name);
mu = d(2:2:end);

inst.name  = name;
inst.Sigma = Sigma;
inst.mu    = mu;
end

function d = scanfile(filename)
fid = fopen(filename);
assert(fid ~= -1, 'cannot open %s', filename);
d = fscanf(fid, '%f');
fclose(fid);
end
