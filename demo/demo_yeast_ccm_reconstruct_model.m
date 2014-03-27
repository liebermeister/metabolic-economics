% DEMO_YEAST_CCM_RECONSTUCT_MODEL - Demo script for reconstruction of enzyme-balanced models

demo_dir = [fileparts(which(mfilename))];

cd(demo_dir)

echo on;
clc
%---------------------------------------------------------------------------------
% DEMO: Economic flux analysis for yeast central metabolism
%
% In this script, we build an enzyme-balanced kinetic model in the following way:
%
% 1.  Determine flux distribution using principle of minimal fluxes
% 2a. Choose thermodynamic forces
% 2b. Choose economic potentials (with the principle of homogeneous costs)
% 3.  Reconstruct an economic-kinetic model realising the flux distribution
%---------------------------------------------------------------------------------
 
% Press key to continue
 
pause
clc
% --------------------------------------------------------------------
% We load a network model of yeast central carbon metabolism.
 
% The file contains the variables: network, network_CoHid, network_CoSplit, v_sign
 
load('/data/yeast_ccm_network.mat');
 
% Press key to continue
 
pause
clc
% --------------------------------------------------------------------
% Phase 1: Determine flux distribution (principle of minimal fluxes)
%          for ATP production, with some predefined flux directions
% ------------------------------------------------------------------
 
% We create structs 'cba_options' and 'cba_constraints' with some default settings:
 
[cba_options, cba_constraints] = cba_default_options(network);
 
% Press key to continue
 
pause
clc
% ------------------------------------------------------------------
% We adjust the settings for flux constraints  
 
cba_constraints.v_sign = v_sign;
 
% We define some reactions to be inactive 
% (strings like 'R00253' are reaction names in our network)
 
cba_constraints.v_fix(label_names({'R00253'},network.actions)) = 0;
cba_constraints.v_fix(label_names({'R00258'},network.actions)) = 0;
cba_constraints.v_fix(label_names({'R00114'},network.actions)) = 0;
cba_constraints.v_fix(label_names({'R00243'},network.actions)) = 0;
cba_constraints.v_fix(label_names({'R00341'},network.actions)) = 0;
cba_constraints.v_fix(label_names({'R00344'},network.actions)) = 0;
cba_constraints.v_fix(label_names({'R00711'},network.actions)) = 0;
 
% We fix one flux direction 
 
cba_constraints.v_sign(label_names({'R00342'},network.actions)) = 1;
 
% We set an upper bound for one reaction 
 
cba_constraints.v_max(label_names({'Oxphos__NADH__irrev__ATP'},network.actions)) = 2;
 
% Press key to continue
 
pause
clc
% ------------------------------------------------------------------
% Next, we adjust the settings for the metabolic objective
 
% Our metabolic objective is the net production of ATP. 
% We declare this by a setting a vector zx with weights 
% for all external production rates
 
zx = zeros(nm,1); 
zx(label_names({'ATP'},network.metabolites)) = 1;
 
% Using this vector, we define the marginal gain vectors (z_ext and zv) 
 
cba_constraints.z_int = 0 * cba_constraints.z_int;
cba_constraints.z_ext = zx(find(network.external));
cba_constraints.zv    = network.N' * zx;
 
% Press key to continue
 
pause
clc
% ------------------------------------------------------------------
% Now we determine the stationary fluxes
 

% 1. Flux Balance Analysis
 
[v_fba,f_benefit] = fba(network,cba_constraints);
 
 
% 2. Fix objective achieved in FBA and run Flux Minimisation (PMF)
 
f_benefit = cba_constraints.zv'*v_fba; 
 
v = pmf(network,cba_constraints,f_benefit,v_fba);
 
 
% 3. Set very small fluxes to 0 and make the flux mode stationary
 
v(abs(v) < 10^-5 *max(abs(v))) = 0;
 
v = project_fluxes(network.N,find(network.external), v,[],sign(v),struct('method','euclidean'));
 
 
% Press key to continue
 
pause
clc
% ---------------------------------------------------------------
% Just to be sure, we apply a function that would remove
% thermodynamically infeasible cycles from our flux mode
 
cba_constraints.ind_ignore = label_names({'Biomass_production'},network.actions);
 
[v,C] = eba_make_feasible(v, network. N, 'loose', nan, cba_constraints.ind_ignore);
 
% Press key to continue
 
pause
clc
% --------------------------------------------------------------------
% Phase 2: Choose chemical and economic potentials
 
% --------------------------------------------------------------------
% We adjust the settings for choosing the chemical potentials
 
cba_constraints.dmu_limit     = 10;
cba_constraints.dmu_limit_min = 2;
cba_constraints.mu_min        = -20 * ones(size(cba_constraints.mu_min));
cba_constraints.mu_max        =  20 * ones(size(cba_constraints.mu_min));
cba_constraints.dmu_min       = -20 * ones(size(cba_constraints.dmu_min));
cba_constraints.dmu_max       =  20 * ones(size(cba_constraints.dmu_min));
cba_constraints.rho           = 100;
  
% Press key to continue
 
pause
clc
% --------------------------------------------------------------------
% We choose the chemical potentials
 
[mu, success_flag] = sample_feasible_mu(network.N,find(network.external),v,cba_constraints,cba_options,'sample',1);
  
% Press key to continue
 
pause
clc
% --------------------------------------------------------------------
% We choose the economic potentials
 
[w, delta_w, y, zx] = cba_homogeneous_cost(network, v, cba_constraints);
 
% Press key to continue
 
pause
clc
% --------------------------------------------------------------------
% Phase 3: Reconstruct an enzyme-balanced model
 
% We set metabolite levels. For the sake of this example, we simply use random values
  
c = 1+5*rand(size(mu));
 
% Now we set some options for the model reconstruction ..
 
cba_options.check_curvatures = 0; 
cba_constraints.Q_ext        = [];
cba_constraints.hu           = ones(size(network.actions)); 
 
% .. and reconstruct the model
 
[network_new, res, cba_constraints_new] = cba_reconstruct_model(network, v, mu, cba_constraints, cba_options, y, w, c);
  
% Press key to continue
 
pause
clc
% ------------------------------------------------------------------
% Finally, we plot some of the reconstructed quantities 
 
% 1. Network and external metabolites
 
figure(1); clf; 
netgraph_concentrations(network_CoSplit,network.external,[],1,struct('actprintnames',1));
  
% Press key to continue
pause
 
% 2. Metabolic fluxes
 
figure(2); clf; 
netgraph_concentrations(network_CoSplit,-network.external,v,1,struct('actstyle','none','arrowsize',0.03));
  
% Press key to continue
pause
 
% 3. Predefined flux signs
 
figure(3); clf; 
netgraph_concentrations(network_CoHid,[],v_sign,1,struct('actstyle','fixed','arrowsize',0.03,'actprintnames',1));
  
% Press key to continue
pause
 
% 4. The local ATP production (contribution of reactions to the metabolic gain)
 
figure(4); clf; 
netgraph_concentrations(network_CoSplit,[],cba_constraints.zv .* v,1,struct('arrowstyle','none'));
  
% Press key to continue
pause
 
% 5. The chemical potentials
 
figure(5); clf;
netgraph_concentrations(network_CoHid,mu,[-network.N'*mu].*[v~=0],1,struct('actstyle','none','arrowsize',0.01));
  
% Press key to continue
pause
 
% 6. The economic potentials
 
figure(5); clf;
netgraph_concentrations(network_CoHid,w,[-network.N'*w].*[v~=0],1,struct('actstyle','none','arrowsize',0.01));
 
% Press key to continue
pause
clc
% That was it - we built a nice kinetic model!
 
% All variables are still in your workspace.
 
% Enjoy working with elasticity sampling!
 
% Press key to finish
pause
return