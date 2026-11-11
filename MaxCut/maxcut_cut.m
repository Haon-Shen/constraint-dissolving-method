function [cut, w] = maxcut_cut(x, L)
% MAXCUT_CUT  Cut obtained from W = reshape(x,N,N) by w = sign(v), v a
% leading eigenvector of W, and its weight 1/4*w'*L*w.
N = size(L, 1);
W = reshape(x, N, N);
[V, d] = eig((W + W')/2, 'vector');
[~, i] = max(d);
w = sign(V(:, i));
w(w == 0) = 1;
cut = w'*L*w/4;
end
