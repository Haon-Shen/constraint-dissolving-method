function [L, A] = maxcut_read(file)
% MAXCUT_READ  Read a graph in rudy format (first line: N, number of edges;
% then one line "i j weight" per edge, nodes numbered from 1) and return the
% Laplacian L = Diag(A*e) - A and the weighted adjacency matrix A.
fid = fopen(file, 'r');
assert(fid >= 3, 'cannot open %s', file);
hdr = fscanf(fid, '%d', 2);
E   = fscanf(fid, '%f', [3, Inf])';
fclose(fid);
N = hdr(1);
assert(size(E, 1) == hdr(2), 'edge count mismatch in %s', file);
A = sparse(E(:,1), E(:,2), E(:,3), N, N);
A = full(A + A');
L = diag(sum(A, 2)) - A;
end
