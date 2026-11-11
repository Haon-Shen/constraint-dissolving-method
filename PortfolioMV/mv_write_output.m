function mv_write_output(results, resdir)
% MV_WRITE_OUTPUT  Plain-text summary and LaTeX table rows for the results
% of run_mv. The figures are drawn by plot_mv_paper.

alg = {'CD', 'ALM', '$1/N$'};
algtex = {'CD', 'ALM~\cite{jia2023augmented}', '$1/N$'};   % labels in the LaTeX tables

% ---------- plain text ----------
fid = fopen(fullfile(resdir, 'mv_summary.txt'), 'w');
for out = [1, fid]
	p = results(1).par;
	fprintf(out, '\nSparse mean-variance problems; CD: c = %g, alpha = %g, tau = %g, eta = %g;  ALM: tau = %g, eta = %g;  eps_out = %g\n', ...
		results(1).c, p.alpha, p.tau, p.eta, p.alm_tau, p.alm_eta, p.eps_out);
	fprintf(out, '\nData sets\n%-12s %5s %5s %5s %6s %4s %8s %9s %12s\n', 'data set', 'n', 'rem.', 'T', 'T_est', 'K', 'windows', 'rho (%)', 'nnz relaxed');
	for r = results
		fprintf(out, '%-12s %5d %5d %5d %6d %4d %8d %9.4f %12d\n', r.name, r.n, r.n_removed, r.T, r.T_est, r.K, r.metrics.comp.windows, r.rho1, r.nnz_relaxed1);
	end
	fprintf(out, '\nIn-sample (first window)\n%-12s %-5s %6s %7s %9s %10s %11s %10s %9s %9s %9s %9s %9s\n', 'data set', 'alg', 'outer', 'inner', 'time(ms)', 'f', 'mu''x-rho', 'Mean(1e-2)', 'Std', 'SR', 'ER', 'RR', 'nnz');
	for r = results
		m = r.metrics;
		for k = 1:2
			o = m.comp.first(k);
			fprintf(out, '%-12s %-5s %6d %7d %9.2f %10.6f %+11.2e %10.4f %9.4f %9.4f %9.4f %9.4f %9d\n', r.name, alg{k}, o.outer, o.inner, 1e3*o.time, o.fval, o.ret_slack, ...
				100*m.in.Mean(k), m.in.Std(k), m.in.SR(k), m.in.ER(k), m.in.RR(k), o.nnz);
		end
	end
	fprintf(out, '\nRolling windows: solver comparison (means over all windows)\n%-12s %8s | %5s %5s %5s | %7s %7s | %8s %8s | %9s %9s | %7s %7s | %9s %9s | %s\n', 'data set', 'windows', ...
		'Cwin', 'Awin', 'tie', 'outC', 'outA', 'innC', 'innA', 'ms C', 'ms A', 'okC', 'okA', 'activeC', 'activeA', 'ties with same support / non-ties with same support');
	for r = results
		q = r.metrics.comp;  tie = abs(q.rel) <= 1e-4;
		fprintf(out, '%-12s %8d | %5d %5d %5d | %7.2f %7.2f | %8.1f %8.1f | %9.2f %9.2f | %7d %7d | %9d %9d | %d / %d   (max |rel| on ties %.1e, rel range [%+.3f, %+.3f])\n', r.name, q.windows, ...
			q.cd_wins, q.alm_wins, q.ties, q.outer, q.inner, 1e3*q.time, q.success, q.active, sum(q.same_support & tie), sum(q.same_support & ~tie), max([0, abs(q.rel(tie))]), min(q.rel), max(q.rel));
	end
	fprintf(out, '\nOut-of-sample\n%-12s %-5s %11s %9s %9s %9s %12s\n', 'data set', 'alg', 'Mean(1e-2)', 'Std', 'SR', 'RR', 'final wealth');
	for r = results
		m = r.metrics;
		for k = 1:3
			fprintf(out, '%-12s %-5s %11.4f %9.4f %9.4f %9.4f %12.4f\n', r.name, strrep(alg{k}, '$', ''), 100*m.out.Mean(k), m.out.Std(k), m.out.SR(k), m.out.RR(k), m.wealth(end,k));
		end
	end
end
fclose(fid);

% ---------- LaTeX tables ----------
% The paper follows the order: data sets, in-sample analysis (first window,
% index data sets), out-of-sample analysis (index data sets), market regimes
% (S&P500 data sets).
reg  = arrayfun(@(r) getfield(mv_dataset_info(r.name), 'regime'), results); %#ok<GFLD>
part = {results(~reg), results(reg)};
fid = fopen(fullfile(resdir, 'mv_tables.tex'), 'w');
fprintf(fid, '%% Table D (data sets): data set & market conditions & stocks & total time period & out-of-sample time period \\\\\n');
for r = results
	info = mv_dataset_info(r.name);
	fprintf(fid, '%s & %s & $%d$ & %s & %s \\\\\n', r.label, info.market, r.n, info.total, info.oos);
end
fprintf(fid, '\n%% Table A (in-sample, first window): data set & Alg. & Outer & Inner & Time (ms) & f & Mean (1e-2) & Std & SR & RR \\\\\n');
for r = part{1}
	m = r.metrics;
	for k = 1:2
		o = m.comp.first(k);
		if k == 1, lead = r.label; else, lead = ''; end
		fprintf(fid, '%s & %s & $%d$ & $%d$ & $%.1f$ & $%.4f$ & $%.3f$ & $%.4f$ & $%.4f$ & $%.2f$ \\\\\n', lead, algtex{k}, ...
			o.outer, o.inner, 1e3*o.time, o.fval, 100*m.in.Mean(k), m.in.Std(k), m.in.SR(k), m.in.RR(k));
	end
end
titles = {'Table B (rolling windows, out-of-sample analysis)', 'Table C (rolling windows, market regimes)'};
for g = 1:2
	fprintf(fid, '\n%% %s: data set & windows & Alg. & Outer & Inner & Time (ms) & N_low & Mean_p (1e-2) & Std_p & SR_p & RR_p \\\\\n', titles{g});
	for r = part{g}
		m = r.metrics;  q = m.comp;  wins = [q.cd_wins, q.alm_wins];
		for k = 1:3
			if k == 1, lead = sprintf('%s & $%d$', r.label, q.windows); else, lead = ' & '; end
			if k <= 2
				mid = sprintf('$%.2f$ & $%.1f$ & $%.1f$ & $%d$', q.outer(k), q.inner(k), 1e3*q.time(k), wins(k));
			else
				mid = '-- & -- & -- & --';
			end
			fprintf(fid, '%s & %s & %s & $%.3f$ & $%.4f$ & $%.4f$ & $%.2f$ \\\\\n', lead, algtex{k}, mid, ...
				100*m.out.Mean(k), m.out.Std(k), m.out.SR(k), m.out.RR(k));
		end
	end
end
fclose(fid);

% ---------- totals quoted in the text ----------
fid = fopen(fullfile(resdir, 'mv_summary.txt'), 'a');
names = {'out-of-sample analysis (index data sets)', 'market regimes (S&P500 data sets)', 'all data sets'};
sets  = {part{1}, part{2}, results};
fprintf(fid, '\nMean time per oracle evaluation (microseconds), CD / ALM / ratio\n');
for r = results
	q = r.metrics.comp;  t = 1e6*q.time./q.nfev;
	fprintf(fid, '%-12s %7.1f %7.1f %6.2f\n', r.name, t, t(1)/t(2));
end
fprintf(fid, '\nTotals\n');
for g = 1:3
	if isempty(sets{g}), continue; end
	q = arrayfun(@(r) r.metrics.comp, sets{g});
	ok = reshape([q.success], 2, []);  ou = reshape([q.outer], 2, []);
	fprintf(fid, '%-42s windows %4d, runs %4d (successful %4d) | CD lower %3d, ALM lower %3d, ties %3d | outer CD %.2f-%.2f, ALM %.2f-%.2f\n', ...
		names{g}, sum([q.windows]), 2*sum([q.windows]), sum(ok(:)), sum([q.cd_wins]), sum([q.alm_wins]), sum([q.ties]), ...
		min(ou(1,:)), max(ou(1,:)), min(ou(2,:)), max(ou(2,:)));
end
fclose(fid);
end
