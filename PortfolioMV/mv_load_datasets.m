function DS = mv_load_datasets(names)
% MV_LOAD_DATASETS  Stock return data sets from the folder data/ (daily
% returns): NASDAQ100 and FTSE100 of Wu et al. (EJOR 2024), file
% init_stock_data.mat, and three S&P500 data sets (declining market,
% COVID-19, rising market), files SnPday_*.mat. names is a cell array of
% the identifiers 'nasdaq100', 'ftse100', 'sp_bear', 'sp_covid', 'sp_bull'
% (default: all).
%
% The data files are not distributed with the code (see data/README.txt).
% init_stock_data.mat contains the matrices NASDAQ100_DAY and FTSE100_DAY,
% and SnPday_Bear.mat, SnPday_COVID_19.mat and SnPday_Bull.mat contain the
% matrices SnPday_Bear, SnPday_COVID_19 and SnPday_Bull. Each matrix holds
% daily simple returns (rows: periods, columns: stocks).
%
% Each entry of DS has the fields
%   name     - identifier
%   label    - name used in tables and figures
%   R         - T x n matrix of returns in percent (rows: periods, columns: stocks)
%   n_removed - number of removed stocks (no price at the beginning of the period)
%   NY        - number of periods per year (250)
%   T_est     - length of the estimation window, round(0.7*T)
%   K         - cardinality bound, min(ceil(n/10), 10)
%
% The first T_est rows form the in-sample data.

datadir = fullfile(fileparts(mfilename('fullpath')), 'data');
files = {'init_stock_data.mat', 'SnPday_Bear.mat', 'SnPday_Bull.mat', 'SnPday_COVID_19.mat'};
missing = files(~cellfun(@(f) isfile(fullfile(datadir, f)), files));
if ~isempty(missing)
	error('mv:nodata', ['Stock return data not found in %s (missing: %s).\n', ...
		'These files are not distributed with the code, see data/README.txt.'], ...
		datadir, strjoin(missing, ', '));
end
A = load(fullfile(datadir, 'init_stock_data.mat'));
B = load(fullfile(datadir, 'SnPday_Bear.mat'));
C = load(fullfile(datadir, 'SnPday_Bull.mat'));
D = load(fullfile(datadir, 'SnPday_COVID_19.mat'));

sets = {
	'nasdaq100', 'NASDAQ100',     A.NASDAQ100_DAY,   250
	'ftse100',   'FTSE100',       A.FTSE100_DAY,     250
	'sp_bear',   'S\&P500-D',     B.SnPday_Bear,     250
	'sp_covid',  'S\&P500-C',     D.SnPday_COVID_19, 250
	'sp_bull',   'S\&P500-R',     C.SnPday_Bull,     250
	};
if nargin < 1, names = sets(:,1)'; end

DS = struct('name', {}, 'label', {}, 'R', {}, 'n_removed', {}, 'NY', {}, 'T_est', {}, 'K', {});
for i = 1:numel(names)
	k = find(strcmp(sets(:,1), names{i}));
	assert(~isempty(k), 'unknown data set %s', names{i});
	DS(i).name  = sets{k,1};
	DS(i).label = sets{k,2};
	R = 100*sets{k,3};
	% Stocks without a price at the beginning of the period are stored with
	% zero returns until their first quotation; they are removed.
	first = zeros(1, size(R, 2));
	for j = 1:size(R, 2)
		q = find(R(:,j) ~= 0, 1);
		if isempty(q), q = size(R, 1) + 1; end
		first(j) = q;
	end
	keep = first <= 5;
	DS(i).R         = R(:, keep);
	DS(i).n_removed = sum(~keep);
	DS(i).NY        = sets{k,4};
	DS(i).T_est     = round(0.7*size(R, 1));
	DS(i).K         = min(ceil(size(DS(i).R, 2)/10), 10);
end
end
