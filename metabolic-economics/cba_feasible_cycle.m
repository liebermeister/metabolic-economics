function [isFeasible, C_infeasible, C_cba, ind_non_orthogonal, n_wasteful] = cba_feasible_cycle(v, N, external, cba_constraints, network, flag_test_wasteful, flag_test_v_sign)

% CBA_FEASIBLE_CYCLE - Check economical flux distribution (non-beneficial cycle criterion)
%
% [isFeasible,C_infeasible,C_cba,ind_non_orthogonal] = cba_feasible_cycle(v, N, external, cba_constraints, network, flag_test_wasteful, flag_test_v_sign)
%
% Calculation:
%  - reduce problem to active subnetwork
%  - check sign-orthogonality on futile cycles


eval(default('flag_test_wasteful','0','flag_test_v_sign','1'));

cba_constraints = cba_update_constraints(cba_constraints,N(find(external),:));
zv              = cba_constraints.zv;

if zv'*v <10^-8, warning('Flux distribution is not beneficial'); end

ind_int = find(external==0); 
N_int   = N(ind_int,:);


% -------------------------------------------------------------

[v_act, N_int_act, Es, nn_act, cba_constraints_act, ind_act] = cba_reduce_to_active_subnetwork(v,N_int,[],network,cba_constraints);

constraint_matrix = cba_constraints_act.zv;
ii =  find(isfinite( cba_constraints_act.v_sign));
nr = size(constraint_matrix,1);

if length(ii),
  display('Treating flux bounds as dual variables like flux benefits');
  for it = 1:length(ii);
    constraint_matrix = [ constraint_matrix, zeros(nr,1)];
    constraint_matrix(ii(it),end) = cba_constraints_act.v_sign(ii(it));
  end
end

C = network_efmtool(nn_act, 'internal', [], constraint_matrix);

ind_non_orthogonal = [];
isFeasible         = 1; 

if size(C), [isFeasible,ind_non_orthogonal] = EBA_orth(sign(v_act),C); end


% -------------------------------------------------------------
% put fluxes back into non-reduced network

C_cba         = zeros(length(v),size(C,2));
C_cba(ind_act,:) = C;

C_infeasible = C_cba(:,ind_non_orthogonal);

if flag_test_v_sign,
  vs = cba_constraints.v_sign;
  if size(C_infeasible,2),
    ind_correct_sign = find(sum(sign(diag(vs(find(isfinite(vs)))) * C_infeasible(find(isfinite(vs)),:) )==-1) ==0);
    C_infeasible = C_infeasible(:,ind_correct_sign);
    ind_non_orthogonal = ind_non_orthogonal(ind_correct_sign);
  else,      ind_correct_sign=[];
     ind_non_orthogonal = [];
  end
end


% -------------------------------------------------------------
% sort infeasible cycles by size

[dum,order]        = sort(sum(abs(sign(C_infeasible)),1));
C_infeasible       = C_infeasible(:,order);
ind_non_orthogonal = ind_non_orthogonal(order);


% test if flux differs from elementary wasteful modes in at least on reaction sign

if flag_test_wasteful,
  C = network_efmtool(nn_act, 'internal', []);
  C = C(:,find(cba_constraints_act.zv'*C<0));
  ind_wasteful = find(sum(diag(sign(v_act))*sign(C)==-1,1)==0);
  n_wasteful = length(ind_wasteful>0); 
end

