function par = default_params()
% DEFAULT_PARAMS  Parameters shared by all experiments. The experiment
% drivers (run_*.m) override the parameters of CD (tau, eta, alpha, beta0)
% and, where stated, eps_out.
%
% The parameters of ALM and of the spectral gradient method are those of
% [MP] = [Jia, Kanzow, Mehlitz, Wachsmuth, Math. Program. 199 (2023),
% Section 6]; those of PD are those of [KL] = [Kanzow, Lapucci, Comput.
% Optim. Appl. 85 (2023)].

% ----- outer loop (all methods) -----
par.eps_out   = 1e-4;    % outer tolerance on the constraint violation (inf-norm)
par.beta_max  = 1e8;     % upper bound on the penalty parameter
par.max_fun   = 1e7;     % budget of oracle (function/gradient) evaluations per run
par.max_outer = Inf;     % optional cap on the number of outer iterations

% ----- penalty update of CD (set for each problem class by the drivers) -----
par.eta       = 0.8;     % required decrease factor of ||G(x)-s||
par.tau       = 2;       % growth factor of the penalty parameter
par.beta0     = [];      % initial penalty parameter; [] = formula of initial_penalty.m

% ----- constraint dissolving mapping of CD -----
par.alpha     = 0.1;     % regularization in M = M_0 + alpha*||H||^2*I

% ----- penalty update of ALM: values of [MP, Section 6] -----
par.alm_eta   = 0.8;     % required decrease factor of V_rho
par.alm_tau   = 10;      % growth factor of the penalty parameter
par.u_min     = -1e20;   % box for the safeguarded multiplier
par.u_max     =  1e20;

% ----- penalty decomposition method PD: [KL, Algorithms 1 and 2] -----
par.pd_tau      = 1.1;   % growth factor of the penalty parameter
par.pd_ls_gamma = 1e-5;  % Armijo parameter gamma of the x-steps
par.pd_ls_beta  = 0.5;   % backtracking factor of the x-steps
par.pd_ls_rtol  = 10*eps; % line search failed if t*|g'*d| <= pd_ls_rtol*max(1,|q|)

% ----- inner solver of CD and ALM: spectral gradient method [MP, Algorithm 3.1] -----
par.eps_inner = 1e-4;    % inner tolerance at outer iteration k is eps_inner/sqrt(k+1)
par.sg_tau    = 2;       % backtracking factor for the inverse stepsize gamma
par.sg_sigma  = 1e-4;    % Armijo parameter
par.sg_m      = 10;      % nonmonotone memory (Grippo et al.)
par.gamma_min = 1e-10;   % bounds for the spectral (BB) inverse stepsize
par.gamma_max = 1e10;
par.gamma_ls_max = 1e20; % line search is declared failed beyond this value
par.max_inner = Inf;     % max. spectral gradient iterations per subproblem

% ----- output -----
par.verbose   = 0;       % 1: print one line per outer iteration
end
