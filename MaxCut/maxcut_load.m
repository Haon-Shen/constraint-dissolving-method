function S = maxcut_load(par)
% MAXCUT_LOAD  The 130 instances of the rudy collection of the BiqMac
% library (folder data/rudy, optimal values in data/rudy_opt.txt):
% S(i).name, S(i).L (Laplacian), S(i).N, S(i).opt (optimal cut value),
% S(i).x0 (solution of the SDP relaxation, see maxcut_sdp_start) and
% S(i).fsdp (optimal value of the relaxation). The starting points are
% computed on the first call (about 5 minutes) and cached in
% results/maxcut_starts_rudy.mat.
if nargin < 1 || isempty(par), par = default_params(); end
here = fileparts(mfilename('fullpath'));
ddir = fullfile(here, 'data');

fid = fopen(fullfile(ddir, 'rudy_opt.txt'), 'r');
fgetl(fid);
C = textscan(fid, '%s %f');
fclose(fid);
names = C{1};  opt = C{2};

cache = fullfile(here, 'results', 'maxcut_starts_rudy.mat');
if exist(cache, 'file')
	load(cache, 'S');
	assert(isequal({S.name}', names), 'cache %s does not match the instance list', cache);
	return
end

S = struct('name', names, 'L', [], 'N', [], 'opt', num2cell(opt), 'x0', [], 'fsdp', [], 'sdp', []);
for i = 1:numel(S)
	S(i).L = maxcut_read(fullfile(ddir, 'rudy', S(i).name));
	S(i).N = size(S(i).L, 1);
	[S(i).x0, o] = maxcut_sdp_start(S(i).L, par);
	S(i).fsdp = o.fsdp;
	S(i).sdp  = rmfield(o, 'inner_flags');
	fprintf('%-20s N = %3d  SDP %.6g (flag %d, %d inner, %.1f s)\n', S(i).name, S(i).N, o.fsdp, o.flag, o.inner, o.time);
end
if ~exist(fileparts(cache), 'dir'), mkdir(fileparts(cache)); end
save(cache, 'S');
end
