function s = portfolio_summary(R)
% PORTFOLIO_SUMMARY  Summary statistics of a batch (output of portfolio_run_batch).
C = [R.cd];  A = [R.alm];
fc = [C.fval];  fa = [A.fval];
s.tol_tie = 1e-4;                        % relative tolerance for a tie
rel = (fc - fa)./max(1, abs(fa));
okc = [C.flag] == 0;  oka = [A.flag] == 0;
rel(~okc & oka) = Inf;  rel(okc & ~oka) = -Inf;   % a run missing the outer criterion loses
s.cd_wins = sum(rel < -s.tol_tie);
s.alm_wins  = sum(rel >  s.tol_tie);
s.ties      = sum(abs(rel) <= s.tol_tie);
s.mean_rel_diff = mean(rel(isfinite(rel)));             % negative: CD lower on average
s.mean_f    = [mean(fc), mean(fa)];
s.outer     = [mean([C.outer]), mean([A.outer])];
s.inner     = [mean([C.inner]), mean([A.inner])];
s.nfev      = [mean([C.nfev]),  mean([A.nfev])];
s.time      = [mean([C.time]),  mean([A.time])];
s.success   = [sum([C.flag] == 0), sum([A.flag] == 0)];
s.inner_fail = [sum([C.inner_fail]), sum([A.inner_fail])];
end
