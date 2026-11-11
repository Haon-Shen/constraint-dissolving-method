function plot_mv_paper(resdir)
% PLOT_MV_PAPER  Figures of the mean-variance problems: cumulative
% out-of-sample wealth of the portfolios computed by CD and ALM and of the equally weighted
% portfolio 1/N. Reads <resdir>/mv.mat (written by run_mv) and writes
%   <resdir>/fig_mv_out.{pdf,png}      NASDAQ100 and FTSE100,
%   <resdir>/fig_mv_regimes.{pdf,png}  the three S&P500 data sets.
% The data sets are split as in mv_write_output. Both figures carry axis
% labels (the captions do not describe the axes) and are 5.125in wide (the
% text width of the paper), so they are included at their natural size.
% Colors as in the other figures of the paper (Okabe-Ito): vermilion for
% CD, blue for ALM, gray for 1/N.

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(resdir), resdir = fullfile(here, 'results'); end
S = load(fullfile(resdir, 'mv.mat'));
results = S.results;
reg = arrayfun(@(r) getfield(mv_dataset_info(r.name), 'regime'), results);

draw(results(~reg), fullfile(resdir, 'fig_mv_out.pdf'));
draw(results(reg), fullfile(resdir, 'fig_mv_regimes.pdf'));
end

function draw(R, outfile)
% One panel per data set in R, legend on top.
if isempty(R), return; end               % no data set of this group was run
col = [213 94 0; 0 114 178; 140 140 140]/255;
lw = [0.8 0.8 0.6];
fs = 7;                                  % \scriptsize of the 10pt paper
h = 1.65;                                % figure height in inches

fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 5.125 h], 'Color', 'w');
tl = tiledlayout(fig, 1, numel(R), 'TileSpacing', 'compact', 'Padding', 'tight');
tl.OuterPosition = [0 0 1 1 - 0.24/h];   % room for the legend on top
for d = 1:numel(R)
	W = [ones(1, 3); R(d).metrics.wealth];
	T = size(W, 1) - 1;
	ax = nexttile(tl);
	hold(ax, 'on');
	for k = 1:3
		plot(ax, 0:T, W(:, k), '-', 'Color', col(k, :), 'LineWidth', lw(k));
	end
	lo = min(W(:));  hi = max(W(:));
	xlim(ax, [0 T]);  ylim(ax, [lo hi] + 0.08*(hi - lo)*[-1 1]);
	set(ax, 'FontSize', fs, 'TickLabelInterpreter', 'latex', 'TickDir', 'in', ...
		'Box', 'on', 'LineWidth', 0.5, 'Layer', 'top');
	ax.YTickLabel = compose('$%g$', ax.YTick);
	ax.XTickLabel = compose('$%g$', ax.XTick);
	% The LaTeX interpreter renders \& as "&amp;", but {\&} correctly.
	title(ax, strrep(R(d).label, '\&', '{\&}'), 'Interpreter', 'latex', 'FontSize', fs, ...
		'FontWeight', 'normal');
	xlabel(ax, 'Out-of-sample period', 'Interpreter', 'latex', 'FontSize', fs);
	if d == 1
		ylabel(ax, 'Cumulative wealth', 'Interpreter', 'latex', 'FontSize', fs);
	end
end
axl = axes(fig, 'Position', [0 1 - 0.15/h 1 0.15/h], 'Visible', 'off');   % legend holder
hold(axl, 'on');
pl = gobjects(1, 3);
for k = 1:3
	pl(k) = plot(axl, NaN, NaN, '-', 'Color', col(k, :), 'LineWidth', lw(k));
end
legend(axl, pl, {'CD', 'ALM', '$1/N$'}, 'Interpreter', 'latex', 'FontSize', fs, ...
	'Orientation', 'horizontal', 'Box', 'on', 'Location', 'north');

exportgraphics(fig, outfile, 'ContentType', 'vector');
exportgraphics(fig, strrep(outfile, '.pdf', '.png'), 'Resolution', 300);
close(fig);
end
