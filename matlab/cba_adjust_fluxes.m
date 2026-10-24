function v_feasible = cba_adjust_fluxes(v, N, external, cba_constraints, network)

% CBA_ADJUST_FLUXES - Correct a flux mode to make it economical
%
% v_feasible = cba_adjust_fluxes(v, N, external, cba_constraints, network)
%
% Adjust given flux distribution v -> economical flux distribution v_feasible
%
% Input 'network' is only necessary if v contains inactive reactions

if cba_constraints.zv'*v <= 0, 
  
  warning('Flux distribution v is not beneficial; adjustment is not possible');
  v_feasible = []; 

else

  epsilon    = 10^-8;

  %% replace approximate zeros by zero
  [this_v,ind_inactive] = v_exact_zeros(v,N,external,epsilon);
  
  ind_int = find(external==0); 
  N_int   = N(ind_int,:);
  [v_act, N_int_act, Es_act, nn_act, cba_constraints_act, ind_act] = cba_reduce_to_active_subnetwork(v, N_int, [], network, cba_constraints);
  CC = network_efmtool(nn_act, 'internal', [], cba_constraints_act.zv);
  CC(abs(CC)<10^-10) = 0;
  C_cba = zeros(length(v),size(CC,2));
  C_cba(ind_act,:) = CC;
  
  % which of the cycles (columns of C_cba) have the same signs as v on their entire support?
  ind_violate = find( sign(this_v)' * C_cba == sum(abs(C_cba)) );
  C_cba       = C_cba(:,ind_violate);
  
  while length(C_cba(:)),
    fprintf('%d conflicts - ',size(C_cba,2));
    clear kappa score

    for it = 1:size(C_cba,2),
      gamma                   = C_cba(:,it);
      [kappa(it),ind_min(it)] = nanmin(this_v(find(gamma))./gamma(find(gamma)));
      score(it)               = kappa(it).^2 * gamma'*gamma;
    end
    
    [score_min,ind]       = min(score);
    [this_v,ind_inactive] = v_exact_zeros(this_v - kappa(ind) * C_cba(:,ind),N,external,epsilon);

    %% keep only cycles whose support is inside the support of the updated v
    C_cba                 = C_cba(:,sum(abs(C_cba(ind_inactive,:)),1) ==0);

    %% keep only cycles that have the same signs as v on their entire support
    ind_violate           = find(sign(this_v)' * C_cba == sum(abs(C_cba)));
    C_cba                 = C_cba(:,ind_violate);
  end
  
  v_feasible = this_v;
  
end
