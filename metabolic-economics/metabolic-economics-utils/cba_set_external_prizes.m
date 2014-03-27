function z_ext = cba_set_external_prizes(network,list);

% z_ext = cba_set_external_prizes(network,list);
%
% construct the vector of external metabolite values
% cell array 'list' contains, in subsequent pairs, 
% a number of metabolite names and the respective values

names  = list(1:2:end-1);
values = cell2mat(list(2:2:end));

ind = label_names(names,network.metabolites);
z = zeros(length(network.metabolites),1);
z(ind) = values;

z_ext = z(find(network.external)); [1 1 1 1 1 1 1 1]';