function [network, res, cba_constraints] = cba_reconstruct_model(network,v,mu,cba_constraints,cba_options,y,w,c)

% CBA_RECONSTRUCT_MODEL - Build model from economical flux mode
%
% [network, res, cba_constraints] = cba_reconstruct_model(network,v,mu,cba_constraints,cba_options,y,w,c)
%
% The function determines all model quantitaties and the rate laws 
% based on an economical flux mode. 
%
% Important: it is assumed that all enzymes in the model are controllable
%
%
% Input
%   network  - Network structure (as in Metabolic Network Toolbox)
%   v        - Precalculated metabolic fluxes
%   w        - Precalculated economic potentials
%   y        - Precalculated enzyme costs
%   mu       - Precalculated chemical potentials
%   c        - Precalculated metabolite levels
%              
%   For the inputs cba_constraints and cba_options, see cba_default_options
%   y and w together must satisfy the economic balance equation
%
% Output
%   network          - Model with rate laws (in field 'kinetics');
%   res              - All results in matlab struct
%   cba_constraints  - Updated constraints data structure


% ----------------------------------------------------------------------
% restrict model to active submodel (variable names ..._act) and initialise variables

% Q_ext is upposed to have length = #external metab  (or 0, if omitted)
if isempty(cba_constraints.Q_ext), 
  Q_ext = nan * ones(sum(network.external),1); 
else, 
  Q_ext = cba_constraints.Q_ext;
end 

if isfield(network,'kinetics'), network = rmfield(network,'kinetics'); end

[nm,nr] = size(network.N);
ind_int = find(network.external==0);
ind_ext = find(network.external==1);

[v_act, N_int_act, Es_act, network_act, cba_constraints_act, ind_active, ind_met_active] = cba_reduce_to_active_subnetwork(v,network.N(ind_int,:),[],network,cba_constraints);

ind_ext_act = find(network_act.external);
ind_int_act = find(network_act.external==0);
mu_act      = mu(ind_met_active);
y_act       = y(ind_active);
w_act       = w(ind_met_active);
c_act       = c(ind_met_active);
N_act       = network_act.N;
W_act       = network_act.regulation_matrix;
% assume hill coefficients = 1 (already implicitly assumed in 'make_structure_matrices')
[Mplus_act, Mminus_act, Wplus_act, Wminus_act, nm_act, nr_act] = make_structure_matrices(N_act,W_act,ind_ext_act);
M_act       = Mplus_act + Mminus_act;
h_act       = ones(nr_act,1); 
n_ext_act   = sum(network_act.external);
A_act       = - N_act'* mu_act;
zeta_act    = exp(h_act .* A_act / RT);
hu_act      = cba_constraints.hu(ind_active);

% note that here, the vector fc refers to all (not only internal) 
% metabolites, but the entries for external metabolites are ignored

if length(cba_constraints_act.zc),
  fc_act = cba_constraints_act.zc(ind_met_active);
else,
  fc_act = zeros(size(c_act));
end

Q_act_predefined = - c_act .* fc_act;
dum              = nan * ones(nm,1);
dum(find(network.external)) = Q_ext;
Q_act_predefined(ind_ext_act) = dum(ind_met_active(ind_ext_act)); 


% ------------------------------------
% find blocks in link matrix

Nint_act = N_act(find(network_act.external ==0),:); 
[L_act, NR_act, ind_ind_met_act] = reduce_N(Nint_act);
L_act_scaled = diag(1./c_act(ind_int_act)) * L_act * diag(c_act(ind_int_act(ind_ind_met_act)));
L_blocks = matrix_find_blocks(L_act);
ind_dep_met_act = setdiff(1:length(ind_int_act),ind_ind_met_act);

[ni,nj] = size(L_act); 
if ni==nj, 
  display(sprintf('\n  The model does not contain conserved moieties')); 
else
  display('  The model contains conserved moieties'); 
end


% -------------------------------------------------------------
% Investment balance equations; precompute necessary information for each metabolite:
%   my_Q       
%   my_beta_M
%   my_alpha_A      
%   my_beta_I     


for it = 1:nm_act,
  %% reactions affected by the metabolite
  ind_rea = find(abs(M_act(:,it))+abs(Wplus_act(:,it))+abs(Wminus_act(:,it)));
  my_M      =       M_act(ind_rea,it);
  my_Mplus  =  Mplus_act(ind_rea,it);
  my_Mminus = Mminus_act(ind_rea,it);
  my_Wplus  =  Wplus_act(ind_rea,it);
  my_Wminus = Wminus_act(ind_rea,it);
  my_y      =      y_act(ind_rea);
  
  % sign constraint for Q // set sign for external metabolites
  Q_min = -inf;  
  Q_max =  inf;  
  
  %% reactions must not be completely switched off or saturated: 0.05 < beta < 0.95
  %% CHECK THE SAME FOR ALLOSTERIC REGLUATION:: FIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  %% ALSO FURTHER BELOW!!!

  x_min = [ Q_min; 0.05 * ones(length(find([my_M; my_Wplus; my_Wminus])),1)];
  x_max = [ Q_max; 0.95 * ones(length(find([my_M; my_Wplus; my_Wminus])),1)];
  
  if network_act.external(it),
    %% external metabolites
    if isfinite(Q_act_predefined(it)),
      my_Q_ext     = Q_act_predefined(it);
      my_Q_ext_std = 0.1*abs(my_Q_ext);
    else
      net_prod  = N_act(it,:) * v_act;
      my_Q_ext = 0;
      if net_prod < 0, my_Q_ext =  1; end
      if net_prod > 0, my_Q_ext = -1; end
      my_Q_ext_std = 1; 
    end
    x_prior_mean = [my_Q_ext;     1/2 * ones(length(find([my_M; my_Wplus; my_Wminus])),1)];
    x_prior_std  = [my_Q_ext_std;       ones(length(find([my_M; my_Wplus; my_Wminus])),1)];
  else
    %% internal metabolites
    my_Q_given = Q_act_predefined(it);
    x_prior_mean = [my_Q_given;           1/2 * ones(length(find([my_M; my_Wplus; my_Wminus])),1)];
    x_prior_std  = [0.00001 * [1 + my_Q_given]; ones(length(find([my_M; my_Wplus; my_Wminus])),1)];
  end
  
  % thermodynamic part of elasticities
  my_E_T = [zeta_act(ind_rea) .* my_Mplus - my_Mminus] ./ [zeta_act(ind_rea)-1];

  my_ind_M      = find(my_M);       my_n_M      = length( my_ind_M);
  my_ind_Wplus  = find(my_Wplus);   my_n_Wplus  = length( my_ind_Wplus);
  my_ind_Wminus = find(my_Wminus);  my_n_Wminus = length( my_ind_Wminus);
  
  Aeq = full([1, ...
              my_y(my_ind_M)'      .*  my_M(my_ind_M)',   ...
            - my_y(my_ind_Wplus)'  .*  my_Wplus(my_ind_Wplus)',   ...
              my_y(my_ind_Wminus)' .* my_Wminus(my_ind_Wminus)']);
  
  Beq = my_y' * my_E_T;
  
  rr{it}.Aeq          = Aeq;
  rr{it}.Beq          = Beq;
  rr{it}.x_min        = x_min;
  rr{it}.x_max        = x_max;
  rr{it}.x_prior_mean = x_prior_mean;
  rr{it}.x_prior_std  = x_prior_std;
  rr{it}.my_n_M       = my_n_M  ;
  rr{it}.my_n_Wplus   = my_n_Wplus  ;
  rr{it}.my_n_Wminus  = my_n_Wminus ;
  
end


% -------------------------------------------------------------
% Initialise more

Q_act       = nan * ones(nm_act,1);
beta_M_act  = zeros(nr_act,nm_act); beta_M_act(find(M_act)) = nan;
alpha_A_act = zeros(nr_act,nm_act); alpha_A_act(find(Wplus_act)) = nan;
beta_I_act  = zeros(nr_act,nm_act);  beta_I_act(find(Wminus_act)) = nan;


% -------------------------------------------------------------
% Solve investment balance equation for external metabolites

opt = optimset('Display','off','Algorithm','interior-point-convex');

for it = 1:n_ext_act,
  ind   = ind_ext_act(it);
  [x_opt,dum,exitflag] = quadprog(diag(rr{ind}.x_prior_std.^-2), -diag(rr{ind}.x_prior_std.^-2)* rr{ind}.x_prior_mean, [],[],rr{ind}.Aeq,rr{ind}.Beq,rr{ind}.x_min,rr{ind}.x_max,[],opt);
  if exitflag ~=1, exitflag 
    error('Error in optimisation'); 
  end 
  
  my_Q               = x_opt(1);                     x_opt = x_opt(2:end);
  my_beta_M_act      = x_opt(1:rr{ind}.my_n_M);      x_opt = x_opt(rr{ind}.my_n_M+1:end);
  my_alpha_A_act     = x_opt(1:rr{ind}.my_n_Wplus);  x_opt = x_opt(rr{ind}.my_n_Wplus+1:end);
  my_beta_I_act      = x_opt(1:rr{ind}.my_n_Wminus);
  
  Q_act(ind,1)                            = my_Q; 
  beta_M_act(find(M_act(:,ind)),ind)      = my_beta_M_act;
  alpha_A_act(find(Wplus_act(:,ind)),ind) = my_alpha_A_act;
  beta_I_act(find(Wminus_act(:,ind)),ind) = my_beta_I_act;
end


% -------------------------------------------------------------
% Solve investment balance equation for internal metabolites
% go through blocks of link matrix and solve equations for each block
% (for instance, blocks with only one entry: independend metabolites
%  on which no other metabolite depends) 

for it = 1:length(L_blocks),
  ind_L_col       = L_blocks{it}.columns;
  ind_L_row       = L_blocks{it}.rows;
  my_L            = L_act(ind_L_row,ind_L_col);
  my_L_scaled     = L_act_scaled(ind_L_row,ind_L_col);

  my_x_prior_std  = [];
  my_x_prior_mean = [];
  my_Aeq          = [];
  my_Beq          = [];
  my_x_min        = [];
  my_x_max        = [];

  for itt = 1:length(ind_L_row)
    it_int = ind_int_act(ind_L_row(itt));
    my_x_min        = [my_x_min;        rr{it_int}.x_min          ];
    my_x_max        = [my_x_max;        rr{it_int}.x_max          ];
    my_x_prior_std  = [my_x_prior_std;  rr{it_int}.x_prior_std;   ];
    my_x_prior_mean = [my_x_prior_mean; rr{it_int}.x_prior_mean   ];
    my_Aeq          = matrix_add_block(my_Aeq,rr{it_int}.Aeq);
    my_Beq          = [my_Beq;          rr{it_int}.Beq            ];    
  end

  my_Aeq = my_L_scaled' * my_Aeq;
  my_Beq = my_L_scaled' * my_Beq;

  [x_opt,fval,exitflag] = quadprog(diag(my_x_prior_std.^-2), -diag(my_x_prior_std.^-2) * my_x_prior_mean, [],[],my_Aeq,my_Beq,my_x_min,my_x_max,[],opt);
  if exitflag ~=1, error('Error in optimisation'); end 
  
  for itt = 1:length(ind_L_row)
    it_int = ind_int_act(ind_L_row(itt));
    my_Q               = x_opt(1);                        x_opt = x_opt(2:end);
    my_beta_M_act      = x_opt(1:rr{it_int}.my_n_M);      x_opt = x_opt(rr{it_int}.my_n_M+1:end);
    my_alpha_A_act     = x_opt(1:rr{it_int}.my_n_Wplus);  x_opt = x_opt(rr{it_int}.my_n_Wplus+1:end);
    my_beta_I_act      = x_opt(1:rr{it_int}.my_n_Wminus); x_opt = x_opt(rr{it_int}.my_n_Wminus+1:end);
  
    Q_act(it_int,1)                               = my_Q; 
    beta_M_act(find(M_act(:,it_int)),it_int)          = my_beta_M_act;
    alpha_A_act(find(Wplus_act(:,it_int)),it_int) = my_alpha_A_act;
    beta_I_act(find(Wminus_act(:,it_int)),it_int) = my_beta_I_act;
  end

end


% --------------------------------------------------------------------------
% update Q values (in active subnetwork) if necessary

Q_mismatch = norm(Q_act(ind_int_act) - Q_act_predefined(ind_int_act));

if Q_mismatch/length(Q_mismatch) < 10^-8,
  display('  Predefined internal economic loads have been realised');
  Q_act(ind_int_act) = Q_act_predefined(ind_int_act);
else,
  display('  Feasible solution requires change of internal economic loads');
  [Q_act(ind_int_act), Q_act_predefined(ind_int_act)]
end

Q_act(abs(Q_act)<10^-8) = 0;


% --------------------------------------
% Convert saturation values to kinetic constants and build kinetics data struture for active submodel

alpha_M_act = alpha_to_betagamma(beta_M_act);
alpha_I_act = alpha_to_betagamma(beta_I_act);

KM_act = alpha_to_k(alpha_M_act,c_act,h_act);
KA_act = alpha_to_k(alpha_A_act,c_act,h_act);
KI_act = alpha_to_k(alpha_I_act,c_act,h_act);

% set u values y/hu, adjust KV values to yield the right flux

network_act.kinetics.type = 'ms';
network_act.kinetics.u    = y_act ./ hu_act;
network_act.kinetics.c    = c_act;
network_act.kinetics.KA   = KA_act;
network_act.kinetics.KI   = KI_act;
network_act.kinetics.KM   = KM_act;
network_act.kinetics.KV   = ones(nr_act,1);
network_act.kinetics.Keq  = exp(N_act' * [log(c_act)-1/RT*mu_act]);
network_act.kinetics.h    = h_act;

vv = network_velocities(c_act, network_act);

if find(vv ==0), error('zero flux encountered'); end
if find(vv.*v_act<0), error('wrong flux direction'); end

network_act.kinetics.KV = v_act ./ vv .* network_act.kinetics.KV;

% check
% v_act_kinetic = network_velocities(c_act, network_act); [v_act_kinetic,v_act]

u_act = network_act.kinetics.u;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check flux benefit response coefficients (should be equal to hu)
% R_act    = basic_control_analysis(network_act,c_act,struct('used',ones(length(v_act),1)));
% Rfvu_act = [ cba_constraints.zv(ind_active)' * R_act.RJ ]';
% Rfcu_act = [ -[Q_act(ind_int_act)./c_act(ind_int_act)]' * R_act.RS(ind_int_act,:) ]';
% ES_act1 = diag(1./v_act) * R_act.epsilon_1 * diag(c_act);
% figure(18); plot(cba_constraints.hu(ind_active), [Rfvu_act + Rfcu_act] .* [v_act~=0],'.');

% check investment condition:
% fc_int_act_updated = -Q_act(ind_int_act)./c_act(ind_int_act);
% E_T_act = diag(1./[zeta_act-1]) * [diag(zeta_act) * Mplus_act - Mminus_act];
% ES_act  = E_T_act - beta_M_act .* M_act + alpha_A_act .* Wplus_act - beta_I_act .* Wminus_act;
% ES_act_unscaled = diag(v_act) * ES_act  * diag(1./c_act); 
% LHS     = - L_act' *  fc_int_act_updated;
% RHS = L_act_scaled' * ES_act(:,ind_int_act)' * y_act
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------------------------------------------------------------------
% Return to entire model with inactive reactions -> compute w and u 
% construct the kinetics for the entire model ...

% set all metabolites from the inactive subnetwork = external

u = 0 * v;    u(ind_active)     = u_act;  
Q = nan * c;  Q(ind_met_active) = Q_act;  

em = sparse(nr,nm);

% for metabolites in inactive part of the network: 
% remove allosteric interactions, but keep them as reactants (KM values)

network.kinetics                 = set_kinetics(network,'ms');
network.kinetics.u               = u;  
network.kinetics.c               = c;
network.kinetics.KA              = em;
network.kinetics.KI              = em;
network.kinetics.KA(ind_active, ind_met_active) = network_act.kinetics.KA ;  
network.kinetics.KI(ind_active, ind_met_active) = network_act.kinetics.KI ;  
network.kinetics.KM(ind_active, ind_met_active) = network_act.kinetics.KM ;
network.kinetics.KV(ind_active)  = network_act.kinetics.KV;
network.kinetics.Keq(ind_active) = network_act.kinetics.Keq;  
network.kinetics.h(ind_active)   = network_act.kinetics.h;  


% ----------------------------------------------------------------------
% set all metabolites from inactive subnetwork external  

ind_met_nonactive = setdiff(1:nm,ind_met_active);
network.external(ind_met_nonactive) = 1;
network.regulation_matrix(:,ind_met_nonactive) = 0;

all_zext = zeros(nm,1);
all_zext(ind_ext) = cba_constraints.z_ext;
cba_constraints.z_ext = all_zext(find(network.external));

% check: v_kinetic = network_velocities(c,network); [v_kinetic,v]


% ----------------------------------------------------------------------
% .. compute the enzyme response coefficients and ..

R    = basic_control_analysis(network,c);
E_sc = diag(1./v) * R.epsilon_1 * diag(c); % scaled elasticities

% % check elasticities:
% % 1. from reconstructed ms rate law
% E_sc_act_1 = E_sc(ind_active,ind_met_active);  
% % 2. from original beta values
% E_T_act    = diag(1./[zeta_act-1]) * [diag(zeta_act) * Mplus_act - Mminus_act];
% E_sc_act_2 = E_T_act - beta_M_act .* M_act + alpha_A_act .* Wplus_act - beta_I_act .* Wminus_act;


% ----------------------------------------------------------------------
% check the balance equations again

% economic balance equation (left hand side, right hand side)
LHS = v.*[network.N'*w + cba_constraints.z_int];
RHS = y;
mismatch = norm(LHS-RHS);
display(sprintf('  Economic balance equation:   mismatch %f',mismatch))
if mismatch>10^-5,  [ LHS, RHS, LHS-RHS]
end

% investment balance equation (left hand side, right hand side)
% first external, then internal metabolites 
% (only metabolites in active subnetwork are considered)

LHS = [Q(ind_met_active(ind_ext_act)); ...
       L_act' * Q(ind_met_active(ind_int_act));];

RHS = [E_sc(ind_active,ind_met_active(ind_ext_act))' * y(ind_active); ...
       L_act_scaled' * E_sc(ind_active,ind_met_active(ind_int_act))' * y(ind_active)];

mismatch = norm(LHS-RHS);
display(sprintf('  Investment balance equation: mismatch %f',mismatch))
if mismatch>10^-5, 
  [ LHS, RHS, LHS-RHS]
end


% ----------------------------------------------------------------------
% .. check the steady state for stability;
% - eigenvalues must not have a positive real part
% - all complex eigenvalues must have a negative real part

M_eigenvalues = eig(full(R.M));
is_stable = [max(real(M_eigenvalues))<=0] * [sum(real(M_eigenvalues(find(imag(M_eigenvalues))))==0)==0];

switch is_stable,
  case 1, display('  o The steady state is stable');
  case 0, warning('  o The steady state is unstable; please sample again!');
end


% ----------------------------------------------------------------------
% .. check if all active enzymes have a positive influence on the benefit 

fc_updated = - Q ./ c;
fc_updated(find(network.external)) = 0;
fc_updated(isnan(fc_updated)) = 0;

Rfvu = [ cba_constraints.zv' * R.RJ ]';
Rfcu = [         fc_updated' * R.RS ]';
Rfu  = Rfvu + Rfcu;

switch min(sign(Rfu(ind_active))),
  case 1,    display('  o All active enzymes have positive marginal benefits');
  otherwise, error('  Active enzyme with non-positive marginal benefit encountered');
end


% --------------------------------------------------------------------------------
% check eigenvalues of fitness curvature matrix
% - for an enzyme-econonic state, the matrix needs to be negative semidefinite
% - to get unique solutions for differential expression prediction, 
%   it needs to be negative definite

if cba_options.check_curvatures,
  %% FIXES NEEDED:
  %% HERE THE FIRST ORDER IS COMPUTED FOR THE SECOND TIME .. MAYBE COMPUTE ONLY ONCE?????
  %% COMPUTING THE TENSOR TAKES LONG .. JUST COMPUTE THE RELEVANT PART!!! 
  [Ec,Eu,parameters,Ecc,Ecu,Euu] = elasticities(network,c,struct('only_enzyme_levels',1));
  [CJ, CS]     = control_coefficients(network.N, Ec,network.external);
  [RSu, RJu, RSuu, RJuu] = response_coefficients(CS, Ec, Eu, Ecc, Ecu, Euu);
  f_u_active = [cba_constraints.zv' * RJu(1:nr,ind_active)]';
  f_uu_active = squeeze(tensor_product(cba_constraints.zv', RJuu(1:nr,ind_active,ind_active)));

  P_orth = eye(length(ind_active)) - 1/[f_u_active'*f_u_active] * [f_u_active*f_u_active'];
  
  %% check eigenvalues in subspace orthogonal on f_u??  IS THAT CORRECT AT ALL??
  [f_uu_active_eigenvectors, f_uu_active_eigenvalues] = eig(full(  P_orth' * f_uu_active * P_orth));
  [f_uu_active_eigenvalues,order] = sort(diag(f_uu_active_eigenvalues));
  f_uu_eigenvectors = f_uu_active_eigenvectors(:,order);
  display(sprintf('  %d directions with positive benefit curvature',sum(f_uu_active_eigenvalues>0)));
  
  figure(1); netgraph_concentrations(network,[],network.kinetics.u,1,struct('arrowstyle','none'));
  title('Enzyme levels');
  show_cc = nan * ones(nr,1); show_cc(ind_active) = f_u_active;
  figure(2); netgraph_concentrations(network,[],network.kinetics.u,1,struct('arrowstyle','none'));
  title('Enzyme benefit gradient');
  show_cc = nan * ones(nr,1); show_cc(ind_active) = f_uu_eigenvectors(:,end);
  figure(3); netgraph_concentrations(network,[],  show_cc,1,struct('arrowstyle','none'));
  title('Curvature eigenvector');
end


% ---------------------------------------------
% output data structure 'res'

u(v==0) = 0;

res.w   = w;
res.u   = u;
res.Q   = Q;
res.fc_updated = fc_updated;
res.kinetics = network.kinetics;
res.R   = R;
res.Rfu = Rfu;
res.Rfvu = Rfvu;
res.Rfcu = Rfcu;

cba_constraints.zc    = fc_updated;
cba_constraints.Q_ext = Q(find(network.external));
