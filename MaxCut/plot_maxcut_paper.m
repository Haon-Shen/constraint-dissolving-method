function plot_maxcut_paper(resdir, outfile)
% PLOT_MAXCUT_PAPER  Bar chart of the cuts on the rudy collection in the
% style of Kanzow-Krueger-Lehmann, Figure 3: for each p = 87, ..., 100 the
% number of instances on which CD, ALM and PD attain at least p percent
% of the optimal cut. Reads <resdir>/maxcut_rudy_{cd,alm,pd}.mat (written
% by run_maxcut; default resdir: results/; PD only if its file is
% complete) and writes <outfile>
% (default <resdir>/fig_maxcut.pdf) and a .png preview. The figure is
% 5.125in wide (the text width of the paper). Colors as in the other
% figures of the paper (Okabe-Ito): vermilion for CD, blue for ALM, and
% bluish green for PD; every bar has a thin black edge.

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(resdir), resdir = fullfile(here, 'results'); end
if nargin < 2 || isempty(outfile), outfile = fullfile(resdir, 'fig_maxcut.pdf'); end
C = load(fullfile(resdir, 'maxcut_rudy_cd.mat'));
A = load(fullfile(resdir, 'maxcut_rudy_alm.mat'));
assert(isequal({C.res.name}, {A.res.name}));
opt = [C.res.opt];
pct = {100*[C.res.cut]./opt, 100*[A.res.cut]./opt};
labels = {'CD', 'ALM'};
fpd = fullfile(resdir, 'maxcut_rudy_pd.mat');
if exist(fpd, 'file')
	D = load(fpd);
	if numel(D.res) == numel(opt)
		assert(isequal({C.res.name}, {D.res.name}));
		pct{end+1} = 100*[D.res.cut]./opt;  labels{end+1} = 'PD';
	end
end

P = 87:100;                              % all three methods attain >= 87% on every instance
cnt = @(v) arrayfun(@(p) sum(v >= p - 1e-12), P);
N = cell2mat(cellfun(cnt, pct, 'UniformOutput', false)')';

col = [213 94 0; 0 114 178; 0 158 115]/255;
fs = 7;                                  % \scriptsize of the 10pt paper

fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 5.125 2.5], 'Color', 'w');
ax = axes(fig);
hold(ax, 'on');
b = bar(ax, P, N, 1, 'grouped', 'EdgeColor', 'k', 'LineWidth', 0.25);
for k = 1:numel(b)
	b(k).FaceColor = col(k, :);
	text(ax, b(k).XEndPoints, b(k).YEndPoints + 2, compose('%d', N(:, k)), ...
		'Rotation', 90, 'FontSize', 5, 'HorizontalAlignment', 'left', ...
		'VerticalAlignment', 'middle', 'Interpreter', 'latex');
end
xlim(ax, [P(1) - 0.6, P(end) + 0.6]);
ylim(ax, [0 numel(opt) + 25]);
set(ax, 'FontSize', fs, 'TickLabelInterpreter', 'latex', 'TickLength', [0 0], ...
	'XTick', P, 'XTickLabel', compose('$%d$', P), ...
	'YTick', 0:20:numel(opt), 'YTickLabel', compose('$%d$', 0:20:numel(opt)), ...
	'YGrid', 'on', 'GridColor', [0.85 0.85 0.85], 'GridAlpha', 1, ...
	'Box', 'on', 'LineWidth', 0.5, 'Layer', 'bottom');
xlabel(ax, 'Percentage of the optimal cut', 'Interpreter', 'latex', 'FontSize', fs);
ylabel(ax, 'Number of instances', 'Interpreter', 'latex', 'FontSize', fs);
legend(ax, b, labels, 'Interpreter', 'latex', 'FontSize', fs, ...
	'Orientation', 'horizontal', 'Box', 'off', 'Location', 'northoutside');

exportgraphics(fig, outfile, 'ContentType', 'vector');
[~, ~, ext] = fileparts(outfile);
exportgraphics(fig, strrep(outfile, ext, '.png'), 'Resolution', 300);
close(fig);
end
