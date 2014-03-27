function [cba_options,cba_constraints] = cba_default_options(network)

% CBA_DEFAULT_OPTIONS - Default settings for directives in 'cba_options' and 'cba_constraints'
%
% [cba_options, cba_constraints] = cba_default_cba_options(network)
%
% Set default values for structures 'cba_options' and 'cba_constraints' 
%
% cba_constraints.v_fix:         vector predetermined fluxes
% cba_constraints.v_min:         vector of lower bounds
% cba_constraints.v_max:         vector of upper bounds
% cba_constraints.v_sign:        predetermined flux signs
% cba_constraints.v_mean:        vector of data values (mean)
% cba_constraints.v_std:         vector of data values (std dev)
% cba_constraints.ext_sign:      sign vector for external metabolite production
%                                
% cba_constraints.mu_fix:        given mu values        
% cba_constraints.mu_min:        lower bounds for mu values  
% cba_constraints.mu_max:        upper bounds for mu values
% cba_constraints.dmu_fix:       given delta mu values        
% cba_constraints.dmu_min:       lower bounds for delta mu values  
% cba_constraints.dmu_max:       upper bounds for delta mu values  
% cba_constraints.dmu_sign:      upper bounds for delta mu values  
% cba_constraints.ind_controllable: which enzymes are controllable?
%
% cba_constraints.fx             cost function (positive) expansion at u=0 (scalar)
% cba_constraints.fxx            cost function (positive) expansion at u=0 (scalar)
% cba_constraints.hu             cost function (positive) gradient at reference state
%
% cba_options.seed               random seed
% cba_options.compute_mu         'mu' or 'delta_mu'
% cba_options.verbose            (Boolean)
% cba_options.test_eba           1
% cba_options.test_cba           1
% cba_options.kinetic_law        'cs'
% cba_options.cba_conditions     'y'
% cba_options.objective          'fba', 'fit'
% cba_options.check_curvatures   default:1

[nm,nr] = size(network.N);

cba_constraints.v_fix       = nan * ones(nr,1);
cba_constraints.v_min       = - ones(nr,1);
cba_constraints.v_max       =   ones(nr,1);
cba_constraints.v_sign      = nan * ones(nr,1);
cba_constraints.v_mean      = nan * ones(nr,1);
cba_constraints.v_std       = nan * ones(nr,1);
cba_constraints.ext_sign    = nan * ones(nm,1);

cba_constraints.mu_fix      = nan * ones(nm,1);
cba_constraints.mu_min      = - ones(nm,1);
cba_constraints.mu_max      =   ones(nm,1);
cba_constraints.dmu_fix     = nan * ones(nr,1);
cba_constraints.dmu_min     = - ones(nr,1);
cba_constraints.dmu_max     =   ones(nr,1);
cba_constraints.dmu_sign    = nan* ones(nr,1);

cba_constraints.y_min       = []; % rate value
cba_constraints.y_max       = []; % rate value

cba_constraints.w_min       = []; % metabolite value
cba_constraints.w_max       = []; % metabolite value

cba_constraints.z_ext       = []; % benefit for production of external metabolites
cba_constraints.z_int       = zeros(nr,1); % direct benefit for fluxes
cba_constraints.zc          = []; % benefit for concentration of internal metabolites

cba_constraints.zx_scaled_min = 0.001;
cba_constraints.u             = nan * ones(nr,1);

cba_constraints.N_tot       = network.N;

% enzyme cost, expansion at u = 0, costs are counted positive
% h(u) = 0 + fx * u + 1/2 * fxx * u^2
cba_constraints.fx          = .1; % first-order enzyme cost
cba_constraints.fxx         = .1; % second-order enzyme cost

% enzyme cost (costs is counted positive), gradient at reference state
cba_constraints.hu          = [];

cba_constraints.ind_controllable = 1:nr; % controllable enzymes

cba_options.seed             = nan;
cba_options.compute_mu       = 'delta_mu';
cba_options.verbose          = 0;
cba_options.test_eba         = 1;
cba_options.test_cba         = 1;
cba_options.kinetic_law      = 'cs';
cba_options.cba_conditions   = 'y';
cba_options.objective        = 'fba'; % 'fit'
cba_options.check_curvatures = 1; 
