function cba_constraints_act = cba_constraints_reduce_to_active(cba_constraints,ind_active,ind_met_active,network,network_act)

% CBA_CONSTRAINTS_REDUCE_TO_ACTIVE - Remove inactive reactions in cba_constraints
%
% cba_constraints_act = cba_constraints_reduce_to_active(cba_constraints,ind_active,ind_met_active,network,network_act)

cba_constraints_act             = cba_constraints;
cba_constraints_act.v_min       = cba_constraints.v_min(ind_active); 
cba_constraints_act.v_max       = cba_constraints.v_max(ind_active); 
cba_constraints_act.v_sign      = cba_constraints.v_sign(ind_active); 
cba_constraints_act.v_mean      = cba_constraints.v_mean(ind_active); 
cba_constraints_act.v_std       = cba_constraints.v_std(ind_active); 
cba_constraints_act.v_fix       = cba_constraints.v_fix(ind_active); 
cba_constraints_act.ext_sign    = cba_constraints_act.ext_sign(ind_met_active);
cba_constraints_act.dmu_min     = cba_constraints.dmu_min(ind_active); 
cba_constraints_act.dmu_max     = cba_constraints.dmu_max(ind_active); 
cba_constraints_act.dmu_fix     = cba_constraints.dmu_fix(ind_active); 
cba_constraints_act.dmu_sign    = cba_constraints.dmu_sign(ind_active); 
cba_constraints_act.mu_min      = cba_constraints.mu_min(ind_met_active); 
cba_constraints_act.mu_max      = cba_constraints.mu_max(ind_met_active); 
cba_constraints_act.mu_fix      = cba_constraints.mu_fix(ind_met_active); 
cba_constraints_act.z_int       = cba_constraints.z_int(ind_active);
cba_constraints_act.u           = cba_constraints.u(ind_active);
cba_constraints_act.ind_ignore  = nan;    %%%%%%%%%%%% FIX!!!
cba_constraints_act.N_tot       = cba_constraints.N_tot(ind_met_active,ind_active);
cba_constraints_act.zv          = [];

if isfield(cba_constraints_act,'hu'),
  cba_constraints_act.hu       = cba_constraints.hu(ind_active);
end

if isfield(cba_constraints_act,'ind_controllable'),
  cba_constraints_act       = rmfield(cba_constraints_act,'ind_controllable');
end

if length(cba_constraints.z_ext),
  nm = length(network.metabolites);
  dummi = zeros(nm,1);
  dummi(find(network.external)) = cba_constraints.z_ext;
  dummi = dummi(ind_met_active);
  new_external = network.external(ind_met_active);
  cba_constraints_act.z_ext = dummi(find(new_external));
end

cba_constraints_act = cba_update_constraints(cba_constraints_act,network_act.N(find(network_act.external),:));
