function p = proj_sparse_nonneg(x, kappa)
% PROJ_SPARSE_NONNEG  One element of the projection of x onto
%
%       D = {x >= 0 : ||x||_0 <= kappa},
%
% obtained by keeping the kappa largest positive entries of x.
p = zeros(size(x));
[v, I] = maxk(x, kappa);
p(I) = max(v, 0);
end
