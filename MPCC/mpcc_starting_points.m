function X0 = mpcc_starting_points(ntries, seed)
% MPCC_STARTING_POINTS  Unprojected starting points near the Rosenbrock
% valley: column i is (w_i, w_i^2 + 0.1*xi_i) with w_i equally spaced in
% [-1,1] and xi_i ~ Unif[0,1]. The solvers project them onto D.
rng(seed);
w  = linspace(-1, 1, ntries);
xi = rand(1, ntries);
X0 = [w; w.^2 + 0.1*xi];
end
