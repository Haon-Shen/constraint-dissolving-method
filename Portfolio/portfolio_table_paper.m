function portfolio_table_paper(resdir)
% PORTFOLIO_TABLE_PAPER  Table of the interval-budget portfolio problems
% with the orl and pard instances merged: one row per pair (n, kappa),
% summarizing the 20 orl instances in <resdir>/orl/portfolio.mat and the 10
% pard instances in <resdir>/pard/portfolio.mat (both written by
% run_portfolio; default resdir: results/). Writes
%   <resdir>/portfolio_table.tex   rows of the LaTeX table
%   <resdir>/portfolio_table.txt   the same summary with further numbers
%                                  quoted in the text (wins on orl and pard
%                                  separately, ratios, active budget bound)

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(resdir), resdir = fullfile(here, 'results'); end
S = load(fullfile(resdir, 'orl', 'portfolio.mat'));
P = load(fullfile(resdir, 'pard', 'portfolio.mat'));

fid = fopen(fullfile(resdir, 'portfolio_table.tex'), 'w');
fprintf(fid, '%% orl + pard; n & kappa & CD: Outer, Inner, time (ms), mean f & ALM: Outer, Inner, time (ms), mean f & N_CD & N_ALM & N_tie \\\\\n');
ftx = fopen(fullfile(resdir, 'portfolio_table.txt'), 'w');
fprintf(ftx, 'orl + pard, merged per (n, kappa); C = CD, A = ALM\n\n');
fprintf(ftx, '   n kappa | wins C/A/tie (orl; pard) | outer C/A | inner C/A (ratio) | time ms C/A (ratio) | mean f C/A | budget active C/A | ties w/ same support\n');
tot = zeros(3, 3);   % rows: orl, pard, all; columns: CD wins, ALM wins, ties
nrun = 0;  nok = 0;  nbud = 0;  nsame = 0;  ntie = 0;
for j = 1:numel(S.results)
	r = S.results(j);  p = P.results(j);
	assert(r.n == p.n && r.kappa == p.kappa);
	R = [r.runs, p.runs];
	s = portfolio_summary(R);
	fprintf(fid, '$%d$ & $%d$ & $%.2f$ & $%.1f$ & $%.1f$ & $%.3f$ & $%.2f$ & $%.1f$ & $%.1f$ & $%.3f$ & $%d$ & $%d$ & $%d$ \\\\\n', ...
		r.n, r.kappa, s.outer(1), s.inner(1), 1e3*s.time(1), s.mean_f(1), s.outer(2), s.inner(2), 1e3*s.time(2), s.mean_f(2), ...
		s.cd_wins, s.alm_wins, s.ties);

	so = portfolio_summary(r.runs);  sp = portfolio_summary(p.runs);
	tot(1, :) = tot(1, :) + [so.cd_wins, so.alm_wins, so.ties];
	tot(2, :) = tot(2, :) + [sp.cd_wins, sp.alm_wins, sp.ties];
	tot(3, :) = tot(3, :) + [s.cd_wins, s.alm_wins, s.ties];

	% outer criterion, active lower budget bound, ties with the same support
	C = [R.cd];  A = [R.alm];
	nrun = nrun + 2*numel(R);
	nok  = nok + sum([C.flag] == 0) + sum([A.flag] == 0);
	bc = sum(abs([C.budget] - 0.8) <= 1e-4);  ba = sum(abs([A.budget] - 0.8) <= 1e-4);
	nbud = nbud + bc + ba;
	rel = ([C.fval] - [A.fval]) ./ max(1, abs([A.fval]));
	same = arrayfun(@(q) isequal(find(q.x_cd), find(q.x_alm)), R);
	tie = abs(rel) <= s.tol_tie;
	ntie = ntie + sum(tie);  nsame = nsame + sum(same & tie);
	fprintf(ftx, '%4d %5d | %2d/%2d/%2d (%2d/%2d/%2d; %2d/%2d/%2d) | %5.2f/%5.2f | %5.1f/%5.1f (%.2f) | %4.1f/%4.1f (%.2f) | %8.3f/%8.3f | %2d/%2d of %d | %d of %d (same support total %d)\n', ...
		r.n, r.kappa, s.cd_wins, s.alm_wins, s.ties, so.cd_wins, so.alm_wins, so.ties, sp.cd_wins, sp.alm_wins, sp.ties, ...
		s.outer, s.inner, s.inner(1)/s.inner(2), 1e3*s.time, s.time(1)/s.time(2), s.mean_f, bc, ba, numel(R), ...
		sum(same & tie), sum(tie), sum(same));
end
fclose(fid);
fprintf(ftx, '\nTotals C/A/tie: orl %d/%d/%d, pard %d/%d/%d, all %d/%d/%d\n', tot');
fprintf(ftx, 'Runs meeting the outer criterion: %d of %d\n', nok, nrun);
fprintf(ftx, 'Runs with e''x within 1e-4 of 0.8: %d of %d\n', nbud, nrun);
fprintf(ftx, 'Ties with the same support: %d of %d\n', nsame, ntie);
fclose(ftx);
type(fullfile(resdir, 'portfolio_table.txt'));
type(fullfile(resdir, 'portfolio_table.tex'));
end
