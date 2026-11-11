function info = mv_dataset_info(name)
% MV_DATASET_INFO  Descriptive information on a data set of mv_load_datasets
% for the data set table of the paper: market condition and the calendar
% periods (mm/dd/yyyy) of all returns and of the out-of-sample returns.
%
% The periods were obtained by matching the return matrices with dated
% prices. The out-of-sample period starts with return number round(0.7*T)+1.
T = {
	'nasdaq100', '--',               '01/02/2020--09/12/2023', '08/03/2022--09/12/2023'
	'ftse100',   '--',               '01/02/2020--09/12/2023', '07/29/2022--09/12/2023'
	'sp_bear',   'Declining market', '03/12/2020--10/12/2022', '01/03/2022--10/12/2022'
	'sp_covid',  'COVID-19',         '12/06/2017--12/31/2020', '02/03/2020--12/31/2020'
	'sp_bull',   'Rising market',    '09/17/2021--09/30/2024', '11/01/2023--09/30/2024'
	};
k = find(strcmp(T(:,1), name));
assert(~isempty(k), 'unknown data set %s', name);
info.market = T{k,2};  info.total = T{k,3};  info.oos = T{k,4};
info.regime = startsWith(name, 'sp_');       % part "market regimes" of the paper
end
