function [v,N_int,Es,network_act,constraints,ind_active,ind_met_active] = cba_reduce_to_active_subnetwork(v,N_int,Es,network,constraints)

% [v,N_int,Es,network_act,constraints,ind_active,ind_met_active] = cba_reduce_to_active_subnetwork(v,N_int,Es,network,constraints)
%
% Reduce cba problem with flux mode v to a cba problem on the active subnetwork

eval(default('Es','[]','network','[]'));

epsilon_v_off = 10^-8;

ind_active = find(abs(v)>epsilon_v_off);
v          = v(ind_active);
N_int      = N_int(:,ind_active);
N_int      = N_int(find(sum(N_int~=0,2)),:);

if ~isempty(Es),
  Es = Es(ind_active,:);
end

if ~isempty(network),
%  network = network_choose(network,1:length(network.metabolites),ind_active);
  [network_act,ind_met_active] = network_choose(network,[],ind_active);
end

if exist('constraints','var'), 
  constraints = cba_constraints_reduce_to_active(constraints,ind_active,ind_met_active,network,network_act);
else
  constraints = [];
end
