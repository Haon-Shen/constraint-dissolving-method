function plot_portfolio_paper(resdir, outfile)
% PLOT_PORTFOLIO_PAPER  Figure of the interval-budget portfolio problems:
% final objective values of ALM and CD on the orl and pard instances, one
% panel per pair (n, kappa). Reads <resdir>/orl/portfolio.mat and
% <resdir>/pard/portfolio.mat, both written by run_portfolio (default
% resdir: results/), and writes <outfile> (default
% <resdir>/fig_portfolio.pdf) and a .png preview.
%
% The objective values on pard are 2-3 times those on orl, so each panel
% consists of two axes side by side: the 20 orl instances (scale on the
% left) and the 10 pard instances (scale on the right). The vertical axes
% do not start at zero. The figure is 5.125in wide (the text width of the
% paper), so it is included at its natural size. Colors as in the other
% figures of the paper (Okabe-Ito): blue for ALM, vermilion for CD.

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(resdir), resdir = fullfile(here, 'results'); end
if nargin < 2 || isempty(outfile), outfile = fullfile(resdir, 'fig_portfolio.pdf'); end
S = load(fullfile(resdir, 'orl', 'portfolio.mat'));
P = load(fullfile(resdir, 'pard', 'portfolio.mat'));

col_alm  = [0 114 178]/255;
col_cd = [213 94 0]/255;
fs = 7;                                  % \scriptsize of the 10pt paper

fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 5.125 4.9], ...
	'Color', 'w');
tl = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'tight');
tl.OuterPosition = [0 0 1 0.955];        % room for the legend on top

for j = 1:numel(S.results)
	assert(S.results(j).n == P.results(j).n && S.results(j).kappa == P.results(j).kappa);
	tj = tiledlayout(tl, 1, 3, 'TileSpacing', 'none', 'Padding', 'none');
	tj.Layout.Tile = j;
	title(tj, sprintf('$n=%d$, $\\kappa=%d$', S.results(j).n, S.results(j).kappa), ...
		'Interpreter', 'latex', 'FontSize', fs);

	% Instance families are named only below the bottom row of panels.
	bottom = j > numel(S.results) - 3;
	ax1 = nexttile(tj, [1 2]);
	draw(ax1, S.results(j).runs, [5.5 15.5], {'orl-05', 'orl-005'}, 'left');
	plot(ax1, [10.5 10.5], ylim(ax1), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
	ax2 = nexttile(tj);
	draw(ax2, P.results(j).runs, 5.5, {'pard'}, 'right');
	if ~bottom, set([ax1 ax2], 'XTickLabel', {}); end

end
axl = axes(fig, 'Position', [0 0.97 1 0.03], 'Visible', 'off');   % legend holder
hold(axl, 'on');
pl = [patch(axl, NaN, NaN, col_alm, 'EdgeColor', 'k', 'LineWidth', 0.25), ...
	patch(axl, NaN, NaN, col_cd, 'EdgeColor', 'k', 'LineWidth', 0.25)];
legend(axl, pl, {'ALM', 'CD'}, 'Interpreter', 'latex', 'FontSize', fs, ...
	'Orientation', 'horizontal', 'Box', 'off', 'Location', 'north');

exportgraphics(fig, outfile, 'ContentType', 'vector');
[~, ~, ext] = fileparts(outfile);
exportgraphics(fig, strrep(outfile, ext, '.png'), 'Resolution', 300);
close(fig);

	function draw(ax, R, xt, xl, side)
	% Grouped bars of the objective values of ALM and CD on the runs R.
		C = [R.cd];  A = [R.alm];
		f = [[A.fval]; [C.fval]]';
		hold(ax, 'on');
		% Two bars of width w side by side at each instance; MATLAB's grouped
		% bars leave too wide a gap between neighboring pairs.
		w = 0.42;  x = (1:size(f, 1))';
		bar(ax, x - w/2, f(:, 1), w, 'FaceColor', col_alm, 'EdgeColor', 'k', 'LineWidth', 0.25);
		bar(ax, x + w/2, f(:, 2), w, 'FaceColor', col_cd, 'EdgeColor', 'k', 'LineWidth', 0.25);
		lo = min(f(:));  hi = max(f(:));
		st = nicestep((hi - lo)/2.5);
		yl = [st*floor((lo - 0.15*(hi - lo))/st), st*ceil((hi + 0.02*(hi - lo))/st)];
		yt = st*(ceil(yl(1)/st):floor(yl(2)/st));
		xlim(ax, [0.5 size(f, 1) + 0.5]);  ylim(ax, yl);
		set(ax, 'FontSize', fs, 'TickLabelInterpreter', 'latex', 'TickLength', [0 0], ...
			'XTick', xt, 'XTickLabel', xl, 'YAxisLocation', side, ...
			'YTick', yt, 'YTickLabel', compose('$%g$', yt), ...
			'YGrid', 'on', 'GridColor', [0.85 0.85 0.85], 'GridAlpha', 1, ...
			'Box', 'on', 'LineWidth', 0.5, 'Layer', 'bottom');
	end
end

function st = nicestep(d)
% Smallest of 1, 2, 5 times a power of ten that is at least d.
p = 10^floor(log10(d));
v = [1 2 5 10];
st = p * v(find(d <= p*v, 1));
end
