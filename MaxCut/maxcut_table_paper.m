function maxcut_table_paper(resdir)
% MAXCUT_TABLE_PAPER  Table of the MaxCut results on the rudy collection:
% for each class of graphs (10 instances each) and for all 130 instances,
% the average numbers of outer and cumulative inner iterations, the
% average running time and the average cut in percent of the optimal
% value, for CD, ALM and PD. Reads <resdir>/maxcut_rudy_{cd,alm,pd}.mat
% (written by run_maxcut; default resdir: results/; PD only if its file
% is complete) and writes
% <resdir>/maxcut_table.tex (table rows, time in seconds) and
% <resdir>/maxcut_table.txt (the same numbers, time in milliseconds,
% together with evaluations and minimum percentages).

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(resdir), resdir = fullfile(here, 'results'); end
C = load(fullfile(resdir, 'maxcut_rudy_cd.mat'));    C = C.res;
A = load(fullfile(resdir, 'maxcut_rudy_alm.mat'));   A = A.res;
assert(isequal({C.name}, {A.name}));
R = {C, A};  labels = {'CD', 'ALM'};
% tex formats of outer, inner, time (s), cut (%)
fmt = {'$%.2f$ & $%.1f$ & $%.3f$ & $%.2f$', '$%.2f$ & $%.1f$ & $%.3f$ & $%.2f$'};
fpd = fullfile(resdir, 'maxcut_rudy_pd.mat');
if exist(fpd, 'file')
	D = load(fpd);  D = D.res;
	if numel(D) == numel(C)
		assert(isequal({C.name}, {D.name}));
		R{end+1} = D;  labels{end+1} = 'PD';
		fmt{end+1} = '$%.1f$ & $%.0f$ & $%.1f$ & $%.2f$';
	end
end
names = {C.name};
cls = regexprep(names, '\.\d+$', '');
[u, ~, g] = unique(cls, 'stable');
groups = [num2cell(1:numel(u)), {0}];     % 0 = all instances

ft = fopen(fullfile(resdir, 'maxcut_table.tex'), 'w');
fx = fopen(fullfile(resdir, 'maxcut_table.txt'), 'w');
p = C(1).par;
fprintf(ft, '%% class');
fprintf(ft, ' & %s: Outer, Inner, time (s), cut (%%)', labels{:});
fprintf(ft, ' \\\\\n');
fprintf(fx, 'MaxCut, rudy collection; CD: tau %g, eta %g, c %g, alpha %g, beta0 %g\n', ...
	p.tau, p.eta, p.c, p.alpha, p.beta0);
fprintf(fx, 'columns per method (%s): outer, inner, f-ev., time (ms), cut (%%), min (%%)\n', strjoin(labels, ', '));
fprintf(fx, '%-9s %3s', 'class', 'n');
fprintf(fx, '%s', repmat(sprintf(' | %8s %9s %10s %10s %7s %7s', 'outer', 'inner', 'f-ev.', 'time(ms)', 'cut(%)', 'min(%)'), 1, numel(R)));
fprintf(fx, '\n');
for k = 1:numel(groups)
	if groups{k} == 0
		i = true(size(names));  lab = 'all';
		fprintf(ft, '\\midrule\n');
	else
		i = (g' == groups{k});  lab = u{groups{k}};
	end
	fprintf(ft, '%s', texname(lab));
	fprintf(fx, '%-9s %3d', lab, sum(i));
	for m = 1:numel(R)
		s = stats(R{m}(i));
		fprintf(ft, [' & ' fmt{m}], s(1), s(2), s(4)/1e3, s(5));
		fprintf(fx, ' | %8.2f %9.1f %10.1f %10.1f %7.2f %7.2f', s);
	end
	fprintf(ft, ' \\\\\n');
	fprintf(fx, '\n');
end
fclose(ft);  fclose(fx);
end

function s = stats(R)
% outer, inner, evaluations, time (ms), average and minimum cut (%)
pc = 100*[R.cut]./[R.opt];
s = [mean([R.outer]), mean([R.inner]), mean([R.nfev]), 1e3*mean([R.time]), mean(pc), min(pc)];
end

function t = texname(lab)
if strcmp(lab, 'all'), t = 'All'; else, t = ['\texttt{' strrep(lab, '_', '\_') '}']; end
end
