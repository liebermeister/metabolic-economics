function [isFeasible, y, dmu] = cba_feasible_lp(v, network, cba_constraints, cba_options)

% CBA_FEASIBLE_LP - Check flux mode for EFA feasibility and choose enzyme costs and thermodynamic forces
%
% [isFeasible, y, dmu] = cba_feasible(v, network, cba_constraints, cba_options)
%
% Test a flux vector v for feasibility (FBA, EBA, CBA) using linear programming
%
% Consider network with internal stoichiometric matrix N_int:
%  - check if flux mode v is enzyme-optimal 
%  - if so, return a possible vector 'y' of rate costs for this state
%    (y = h_u * u / v), same signs as the rates
%
% Structure 'cba_constraints' must contains the fields:
%
% either    zv   : benefits for fluxes 
% or        z_ext: benefits of network.external metabolites 
%
% N_tot     : entire stoichiometric matrix; required for EBA test
%
% cba_options.test_eba : flag: should thermodynamical feasibility be tested as well?


% ------------------------------------------------------------------
% if necessary, reduce to active subnetwork and call the function

% set all small fluxes to zero
v(abs(v)/max(abs(v))<10^-5)=0;

if sum(v==0),
  N_int = network.N(find(network.external==0),:);
  Es = -network.N';
  [v_act,N_int_act,Es_act,network_act,cba_constraints_act,ind_act,ind_met_act] = cba_reduce_to_active_subnetwork(v,N_int,Es,network,cba_constraints);
  [isFeasible, y_act, dmu_act] = cba_feasible_lp(v_act,network_act,cba_constraints_act,cba_options);
  [nm,nr]             = size(network.N);
  y           = nan *ones(nr,1);
  y(ind_act)  = y_act;  
  dmu            = nan *ones(nr,1);
  if isFeasible,
    dmu(ind_act)   = dmu_act;
  end
  return
end


% ------------------------------------------------------------------
%% initialisation

cba_constraints = cba_update_constraints(cba_constraints,network.N(find(network.external),:));

ind_int = find(network.external==0);
N_int   = network.N(ind_int,:);
isFeasible  = 1;
y         = [];
dmu  = [];
zv        = cba_constraints.zv;
  
if ~isfield(cba_options,'test_eba'), cba_options.test_eba = 1; end

%% test if mode is stationary
  
epsilon_balance = 10^-4;
  
if max(abs(N_int*v))>epsilon_balance,
  isFeasible = 0;  display('Flux mode is not stationary.');
end

% ------------------------------------------------------------------
% first condition: mode yields positive benefit?

v_w     = zv'*v;
  
switch sign(v_w),
  case  0,  isFeasible = 0; display('No benefit produced.');
  case -1,  isFeasible = 0; display('Negative  benefit produced.');
end

% ------------------------------------------------------------------
% second condition: existence of enzyme cost vector y

if isFeasible,  
    
  %% reduce model to active reactions and remove all
  %% irrelevant internal metabolites
  
  nn = network_construct(N_int);
  cc = cba_constraints;
  cc.z_ext = [];
  
  %% test if sign pattern is feasible: search for a vector y such that
  %%
  %%      K' * y = K' * zv
  %% diag(v) * y > 0
  %%
  %% solve this by linear programming:
    
  K  = analyse_N(N_int);    
  nr = size(N_int,2);
    
  epsilon = 10^-4; % threshold enzyme cost
  c       = ones(nr,1);
  h       = epsilon * ones(nr,1);
  G       = diag(v);
  A       = K';
  b       = K' * cc.zv;

  %%y_act = lp236a(-c,-G,-h,A,b);

  if exist('cplexlp','file'),
    opt = optimset('Display','off');
    [y, fval, exitflag] = linprog(-c, -G, -h, A, b,[],[], [], opt);
  else
    opt = cplexoptimset('Display','off');
    [y, fval, exitflag] = cplexlp(-c, -G, -h, A, b,[],[], [], opt);
  end

  if exitflag ~= 1,
    exitflag
    isFeasible = 0;  display('Flux distribution seems to be economically infeasible.');
  end
  
  if cba_options.test_eba,
    [eba_feas,dmu] = eba_feasible_lp(v,cba_constraints.N_tot);
    isFeasible = isFeasible * eba_feas;
  end

end
