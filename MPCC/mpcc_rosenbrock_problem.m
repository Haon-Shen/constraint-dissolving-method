function prob = mpcc_rosenbrock_problem(c)
% MPCC_ROSENBROCK_PROBLEM  Interval-constrained Rosenbrock problem
%
%       min  (1-x1)^2 + 100*(x2-x1^2)^2
%       s.t. 1/2 <= x1 + x2 <= 3/2,   x1 >= 0,  x2 >= 0,  x1*x2 = 0,
%
% in the form  min f(x)  s.t.  G(x) in C,  x in D  with
%
%       G(x) = x1 + x2,   C = [1/2, 3/2],   D = {x >= 0 : x1*x2 = 0},
%
% and the projective mappings
%
%       Q_D(x) = c*diag(x1^2, x2^2),      Q_C(s) = (s - 1/2)*(3/2 - s).
%
% Global minimizer x* = (1/2,0), f = 6.5; local minimizer (0,1/2), f = 26.

if nargin < 1, c = 1; end

lo = 0.5;  up = 1.5;                 % C = [lo, up]

prob.name = 'mpcc_rosenbrock';
prob.n = 2;
prob.m = 1;
prob.c = c;

prob.f     = @(x) (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2;
prob.gradf = @(x) [-2*(1 - x(1)) - 400*x(1)*(x(2) - x(1)^2);
                   200*(x(2) - x(1)^2)];

prob.G     = @(x) x(1) + x(2);
prob.JG    = @(x) [1, 1];
prob.hessG = [];                     % G is affine

prob.projD = @(x) proj_complementarity(x);
prob.projC = @(s) proj_box(s, lo, up);

prob.QD  = @(x, U) c*(x.^2).*U;
prob.dQD = @(x, u1, u2) 2*c*x.*u1.*u2;
prob.QC  = @(s, W) ((s - lo).*(up - s)).*W;
prob.dQC = @(s, w1, w2) (lo + up - 2*s).*w1.*w2;

% reference solutions, used to classify the computed points
prob.x_global = [0.5; 0];
prob.x_local  = [0; 0.5];
end
