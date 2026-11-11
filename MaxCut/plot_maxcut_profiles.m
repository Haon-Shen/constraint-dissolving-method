function plot_maxcut_profiles(resdir, outfile)
% PLOT_MAXCUT_PROFILES  Cumulative curves of the numerical effort on the
% rudy collection in the style of Kanzow-Krueger-Lehmann, Figure 4: for
% each method, the number of instances solved within a given running time,
% number of outer iterations and number of cumulative inner iterations
% (three panels in one row). Reads <resdir>/maxcut_rudy_{cd,alm}.mat
% (written by run_maxcut; default resdir: results/) and writes <outfile> (default
% <resdir>/fig_maxcut_profiles.pdf) and a .png preview. An instance counts
% as solved only if the run meets its outer stopping criterion. The figure
% is 5.125in wide (the text width of the paper) and each panel has the
% aspect ratio 8:5; the figure height 1.5in is chosen so that the panels
% fill the width. Colors as in the other
% figures of the paper (Okabe-Ito): vermilion for CD, blue for ALM.

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(resdir), resdir = fullfile(here, 'results'); end
if nargin < 2 || isempty(outfile), outfile = fullfile(resdir, 'fig_maxcut_profiles.pdf'); end
C = load(fullfile(resdir, 'maxcut_rudy_cd.mat'));  C = C.res;
A = load(fullfile(resdir, 'maxcut_rudy_alm.mat'));   A = A.res;
assert(isequal({C.name}, {A.name}));
n = numel(C);

fields = {'time', 'outer', 'inner'};
labels = {'(a) Running time (s)', '(b) Outer iterations', '(c) Inner iterations'};
col = [213 94 0; 0 114 178]/255;
fs = 7;                                  % \scriptsize of the 10pt paper

fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 5.125 1.5], 'Color', 'w');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'tight');
tl.OuterPosition = [0 0 1 0.9];         % room for the legend on top
for k = 1:numel(fields)
	ax = nexttile(tl);
	hold(ax, 'on');
	h = gobjects(1, 2);
	R = {C, A};
	for m = 1:2
		v = [R{m}.(fields{k})];
		v = sort(v([R{m}.flag] == 0));
		h(m) = stairs(ax, [0, v, max(v)*1.05], [0, 1:numel(v), numel(v)], '-', ...
			'Color', col(m, :), 'LineWidth', 0.8);
	end
	vmax = max([[C.(fields{k})], [A.(fields{k})]]);
	xlim(ax, [0, 1.05*vmax]);  ylim(ax, [0, n + 5]);
	pbaspect(ax, [8 5 1]);                  % each panel 8:5 (width:height)
	set(ax, 'FontSize', fs, 'TickLabelInterpreter', 'latex', ...
		'YTick', 0:20:n, 'YTickLabel', compose('$%d$', 0:20:n), ...
		'XGrid', 'on', 'YGrid', 'on', 'GridColor', [0.85 0.85 0.85], 'GridAlpha', 1, ...
		'Box', 'on', 'LineWidth', 0.5, 'Layer', 'bottom');
	xt = get(ax, 'XTick');
	set(ax, 'XTickLabel', compose('$%g$', xt));
	xlabel(ax, labels{k}, 'Interpreter', 'latex', 'FontSize', fs);
	if k == 1
		ylabel(ax, 'Number of instances', 'Interpreter', 'latex', 'FontSize', fs);
	end
end
axl = axes(fig, 'Position', [0 0.92 1 0.08], 'Visible', 'off');   % legend holder
hold(axl, 'on');
pl = [plot(axl, NaN, NaN, '-', 'Color', col(1, :), 'LineWidth', 0.8), ...
	plot(axl, NaN, NaN, '-', 'Color', col(2, :), 'LineWidth', 0.8)];
legend(axl, pl, {'CD', 'ALM'}, 'Interpreter', 'latex', 'FontSize', fs, ...
	'Orientation', 'horizontal', 'Box', 'off', 'Location', 'north');

exportgraphics(fig, outfile, 'ContentType', 'vector');
[~, ~, ext] = fileparts(outfile);
exportgraphics(fig, strrep(outfile, ext, '.png'), 'Resolution', 300);
close(fig);
end
