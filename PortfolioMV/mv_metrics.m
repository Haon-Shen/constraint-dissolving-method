function M = mv_metrics(ds, W)
% MV_METRICS  Evaluation of the rolling-window results W (output of
% mv_rolling for all windows) of the data set ds, following the criteria of
% Wu et al. (EJOR 2024). Returns are converted from percent to plain numbers.
%
% M.in   - in-sample criteria of the first window, rows: CD, ALM
%          fields Mean, Std, SR, ER, RR
% M.out  - out-of-sample criteria, rows: CD, ALM, 1/N
%          fields Mean, Std, SR, RR and the cumulative wealth (columns of M.wealth)
% M.comp - comparison of the two solvers over all windows

Rf = 0.036;                              % annual risk-free rate
n  = size(ds.R, 2);
yn = ones(n, 1)/n;                       % 1/N portfolio

% ----- in-sample, first window -----
Rest  = ds.R(1:ds.T_est, :)/100;
mu    = mean(Rest, 1)';  Sigma = cov(Rest);
rho   = W(1).rho/100;
X     = [normalize_w(W(1).x_cd), normalize_w(W(1).x_alm)];
for k = 1:2
	y = X(:,k);
	M.in.Mean(k,1) = mu'*y;
	M.in.Std(k,1)  = sqrt(y'*Sigma*y);
	M.in.SR(k,1)   = (M.in.Mean(k) - Rf/ds.NY)/M.in.Std(k);
	M.in.ER(k,1)   = (M.in.Mean(k) - rho)/rho;
	M.in.RR(k,1)   = (yn'*Sigma*yn)/(y'*Sigma*y);
end

% ----- out-of-sample, all windows -----
Kw    = numel(W);
Rout  = vertcat(W.r_out)/100;            % Kw x n
Sout  = cov(Rout);
Yc    = cell2mat(arrayfun(@(w) normalize_w(w.x_cd),   W, 'UniformOutput', false));
Ya    = cell2mat(arrayfun(@(w) normalize_w(w.x_alm),  W, 'UniformOutput', false));
Yn    = repmat(yn, 1, Kw);
Ys    = {Yc, Ya, Yn};
vn    = yn'*Sout*yn;
for k = 1:3
	r = sum(Rout'.*Ys{k}, 1)';           % realized portfolio returns
	M.out.Mean(k,1) = mean(r);
	M.out.Std(k,1)  = std(r);
	M.out.SR(k,1)   = (mean(r) - Rf/ds.NY)/std(r);
	M.out.RR(k,1)   = mean(vn./sum(Ys{k}.*(Sout*Ys{k}), 1));
	M.wealth(:,k)   = cumprod(1 + r);
end

% ----- solver comparison over all windows -----
C = [W.cd];  A = [W.alm];
rel = ([C.fval] - [A.fval])./max(1, abs([A.fval]));
tol = 1e-4;
M.comp.windows   = Kw;
M.comp.cd_wins = sum(rel < -tol);
M.comp.alm_wins  = sum(rel > tol);
M.comp.ties      = sum(abs(rel) <= tol);
M.comp.rel       = rel;
M.comp.same_support = arrayfun(@(w) isequal(find(w.x_cd), find(w.x_alm)), W);
M.comp.outer     = [mean([C.outer]), mean([A.outer])];
M.comp.inner     = [mean([C.inner]), mean([A.inner])];
M.comp.nfev      = [mean([C.nfev]),  mean([A.nfev])];
M.comp.time      = [mean([C.time]),  mean([A.time])];
M.comp.mean_f    = [mean([C.fval]),  mean([A.fval])];
M.comp.success   = [sum([C.flag] == 0), sum([A.flag] == 0)];
M.comp.inner_fail = [sum([C.inner_fail]), sum([A.inner_fail])];
M.comp.active    = [sum(abs([C.ret_slack]) <= 1e-4), sum(abs([A.ret_slack]) <= 1e-4)];
M.comp.first     = [C(1), A(1)];
end

function y = normalize_w(x)
y = x/sum(x);
end
