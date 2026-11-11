function x = proj_rank1_psd(x, N)
% PROJ_RANK1_PSD  A projection of W = reshape(x,N,N) onto
% {W psd : rank(W) <= 1}: max(lambda,0)*v*v' with (lambda,v) a largest
% eigenpair of the symmetric part of W [MP 2023, Proposition 5.7].
W = reshape(x, N, N);
W = (W + W')/2;
[V, d] = eig(W, 'vector');
[lam, i] = max(d);
v = V(:, i);
x = reshape(max(lam, 0)*(v*v'), [], 1);
end
