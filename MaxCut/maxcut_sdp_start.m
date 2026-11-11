function [x0, out] = maxcut_sdp_start(L, par)
% MAXCUT_SDP_START  Starting point: a solution of the SDP relaxation
%
%       min -1/4*trace(L*W)   s.t.  diag(W) = e,  W psd,
%
% computed, as in [MP 2023, Section 6.3], by the safeguarded ALM with
% D = {W psd} (no toolbox needed), starting from the identity matrix.
% L is scaled by 1/max|L_ij|, which does not change the solution set.
% out.fsdp is the
% optimal value -1/4*trace(L*W) of the relaxation (an upper bound on the
% maximum cut up to the tolerance).
if nargin < 2, par = default_params(); end
N    = size(L, 1);
sc   = max(abs(L(:)));
prob = maxcut_problem(L/sc);
prob.projD = @(x) proj_psd(x, N);
[x0, ~, out] = safeguarded_alm(prob, reshape(eye(N), [], 1), par);
out.fsdp = -sc*out.fval;
end

function x = proj_psd(x, N)
W = reshape(x, N, N);
[V, d] = eig((W + W')/2, 'vector');
x = reshape(V*diag(max(d, 0))*V', [], 1);
end
