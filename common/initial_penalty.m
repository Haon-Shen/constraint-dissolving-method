function beta0 = initial_penalty(fx0, viol0)
% INITIAL_PENALTY  Initial penalty parameter of ALM and PD, and of CD
% unless par.beta0 is set [MP, Section 6],
%
%   beta0 = P_[1e-3,1e3]( 10*max{1,f(x0)} / max{1, 0.5*viol0^2} ),
%
% where viol0 = ||G(x0)-s0||_2 with s0 = P_C(G(x0)), i.e. dist_C(G(x0)).
beta0 = proj_box(10*max(1, fx0)/max(1, 0.5*viol0^2), 1e-3, 1e3);
end
