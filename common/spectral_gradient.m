function [z, out] = spectral_gradient(phi, proj, z, tol, par, max_fun)
% SPECTRAL_GRADIENT  Nonmonotone spectral projected gradient method for
%
%       min phi(z)   s.t.  z in K      (K closed, possibly nonconvex),
%
% see [Jia, Kanzow, Mehlitz, Wachsmuth, Math. Program. 2023, Algorithm 3.1].
%
% Input:
%   phi     - handle, [val, grad] = phi(z)
%   proj    - handle, one element of the projection of its argument onto K
%   z       - starting point, must belong to K
%   tol     - tolerance of the termination criterion
%               || gamma*(z_j - z_trial) + grad(z_trial) - grad(z_j) ||_inf <= tol,
%             which is tested at every projected trial point (Line 5 of Alg. 3.1)
%   par     - parameters, see default_params
%   max_fun - (optional) budget of evaluations of phi
%
% Output:
%   z   - computed point (in K)
%   out - struct with fields
%           iter : number of iterations
%           nfev : number of evaluations of phi
%           res  : value of the termination measure at the last trial point
%           flag : 0 tolerance reached, 1 max_inner reached,
%                  2 evaluation budget exhausted, 3 line search failed

if nargin < 6, max_fun = Inf; end

[fz, gz] = phi(z);
nfev  = 1;
iter  = 0;
res   = Inf;
flag  = 0;
fhist = fz;                       % last (at most m+1) accepted function values

while true
	if iter >= par.max_inner, flag = 1; break; end

	% spectral (Barzilai-Borwein) choice of the initial inverse stepsize
	if iter == 0
		gamma = norm(gz, inf);
	else
		gamma = (dz'*dg)/(dz'*dz);
	end
	gamma = proj_box(gamma, par.gamma_min, par.gamma_max);
	fmax  = max(fhist);

	% backtracking with the nonmonotone Armijo rule
	accepted = false;
	while true
		if nfev >= max_fun,           flag = 2; break; end
		if gamma > par.gamma_ls_max,  flag = 3; break; end

		znew = proj(z - gz/gamma);
		[fnew, gnew] = phi(znew);
		nfev = nfev + 1;
		d    = znew - z;

		res = norm(gnew - gz - gamma*d, inf);
		if res <= tol || fnew <= fmax + par.sg_sigma*(gz'*d)
			accepted = true;
			break
		end
		gamma = par.sg_tau*gamma;
	end
	if ~accepted, break; end

	dz = d;  dg = gnew - gz;
	z  = znew;  fz = fnew;  gz = gnew;
	fhist = [fhist(max(1, end-par.sg_m+1):end); fz];
	iter  = iter + 1;

	if res <= tol, break; end
end

out.iter = iter;
out.nfev = nfev;
out.res  = res;
out.flag = flag;
end
