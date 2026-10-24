function result = cba_economic_state(N, external, cba_constraints, cba_options)

% CBA_ECONOMIC_STATE - Find fluxes, chemical and economic potentials satisfying FBA/EBA/CBA constraints
%
% result = cba_economic_state(N, external, cba_constraints, cba_options)
%
% Determine feasible flux vector v, chemical potential differences dmu, 
% and rate values y under EBA and CBA constraints with either
%  - linear FBA objective (zv' * v) 
%  - sum of squared residuals objective sum([[v-v_mean]/v_std].^2)
% plus a small Euclidean norm term for regularisation
%
% Output:
%   result.v      - Flux mode
%   result.dmu    - Chemical potential differences
%   result.y      - (only exists if cba_options.cba_conditions == 'y')
%   result.w      - (only exists if cba_options.cba_conditions == 'w')
%   result.value  - Objective function


% The constraints given by cba_constraints.mu_fix, cba_constraints.dmu_fix are ignored


%  ------------------------------------------------------------
% initialisation

[nm,nr]      = size(N);

epsilon = 10^-3;

cba_constraints.v_min(find(cba_constraints.v_sign==1))  = epsilon;
cba_constraints.v_max(find(cba_constraints.v_sign==-1)) = -epsilon;

if isempty(cba_constraints.dmu_min),  cba_constraints.dmu_min = -10*ones(nr,1); end
if isempty(cba_constraints.dmu_max),  cba_constraints.dmu_max =  10*ones(nr,1); end
if isempty(cba_constraints.y_min),    cba_constraints.y_min   = -10*ones(nr,1); end
if isempty(cba_constraints.y_max),    cba_constraints.y_max   =  10*ones(nr,1); end
if isempty(cba_constraints.w_min),    cba_constraints.w_min   = -10*ones(nm,1); end
if isempty(cba_constraints.w_max),    cba_constraints.w_max   =  10*ones(nm,1); end

zv           = cba_constraints.zv;
ind_v_fix    = find(isfinite(cba_constraints.v_fix));
ind_extsigns = find(isfinite(cba_constraints.ext_sign));
n_int        = sum(external==0); 
N_int        = N(find(external==0),:);
Ktot         = analyse_N(N);
Kint         = analyse_N(N(find(external==0),:));
nKint        = size(Kint,2);
nKtot        = size(Ktot,2);
my_eye       = eye(nr);
regularisation_weight  = 0; % weights for regularisation term


% ------------------------------------------------------------
% vector x contains all (v, dmu, y) or (v, dmu, w)
% indices within vector (either ind_y or ind_w is unnecessary)

ind_v   =           1:nr;
ind_dmu =     nr + (1:nr);
ind_y   = 2 * nr + (1:nr);
ind_w   = 2 * nr + (1:nm);


% ------------------------------------------------------
% cba_constraints:
%
% (where v = x(ind_v), dmu = x(ind_dmu), y = x(ind_y), w = x(ind_w) ) 
%
% limits           : cba_constraints.v_min   <=  x(ind_v)   <= cba_constraints.v_max
%                  : cba_constraints.dmu_min <=  x(ind_dmu) <= cba_constraints.dmu_max
%                  : cba_constraints.y_min   <=  x(ind_y)   <= cba_constraints.y_max
% regularisation   : sum(x.^2) = min !
%
% fixed fluxes     : x(ind_v_fix)          = cba_constraints.v_fix(ind_v_fix)
% stationary fluxes: N_int * x(ind_v )   = 0
% positive objective: zv' * x(ind_v )   > 0
%
% eba, wegscheider : Ktot' * x(ind_dmu)  = 0
% eba signs        : diag(x(ind_v)) * x(ind_dmu) <= - epsilon < 0   where v == 0
%
% variant 1: use rate values y
%   y balance      : Kint' * x(ind_y)    = Kint' * zv
%   y signs        : diag(x(ind_v)) * x(ind_y) >= epsilon > 0   where v == 0
%
% variant 2: use metabolite values w
%   w signs        : diag(x(ind_v)) * [ N' * x(ind_w) + zv] >= epsilon > 0 if v == 0

% ------------------------------------------------------
% rewrite this in the form 
% 
% Aeq * x  = beq
%   A * x <= b
%       C  = cbaf2(x,N)   (for signs)
%     min != cbaf1(x)     (regularisation)


switch cba_options.cba_conditions, 
  
  case 'y', 

    display('Optimisation with rate values y');

    %% Equality cba_constraints:
    %%  v_fix
    %%  stationarity
    %%  wegscheider for dmu
    %%  balance for y
    
    Aeq  = double([my_eye(ind_v_fix,:), zeros(length(ind_v_fix),2*nr)      ; ...
                   N_int,               zeros(n_int,2*nr)                  ; ...
                   zeros(nKtot,nr),     Ktot',            zeros(nKtot,nr)  ; ...
                   zeros(nKint,2*nr),                     Kint'             ]);

    beq  = double([cba_constraints.v_fix(ind_v_fix);... 
                   zeros(n_int,1)              ;...
                   zeros(nKtot,1)              ;...
                   Kint' * zv ]);

    xmin   = [cba_constraints.v_min; cba_constraints.dmu_min; cba_constraints.y_min];
    xmax   = [cba_constraints.v_max; cba_constraints.dmu_max; cba_constraints.y_max];
    zvfull = [cba_constraints.zv ; zeros(nr,1); zeros(nr,1)];

    %% inequality cba_constraints: 
    %% positive objective
    %% external signs
    
    A = [-zv', zeros(1,2*nr) ; ...
         double([-diag(cba_constraints.ext_sign(ind_extsigns)) * N(ind_extsigns,:), ...
                 zeros(length(ind_extsigns),2*nr)])];
    
    b = [0; zeros(length(ind_extsigns),1)];
    
  case 'w',  

    display('Optimisation with metabolite values w');

    %% Equality cba_constraints:
    %%  v_fix
    %%  stationarity
    %%  wegscheider for dmu
    
    Aeq  = double([my_eye(ind_v_fix,:),    zeros(length(ind_v_fix),nr+nm) ; ...
                   N_int,                 zeros(n_int,nr+nm)              ; ...
                   zeros(nKtot,nr),       Ktot',          zeros(nKtot,nm)]);
    
    beq  = double([cba_constraints.v_fix(ind_v_fix);... 
                   zeros(n_int,1)           ;...
                   zeros(nKtot,1)          ]);

    xmin   = [cba_constraints.v_min; cba_constraints.dmu_min; cba_constraints.w_min];
    xmax   = [cba_constraints.v_max; cba_constraints.dmu_max; cba_constraints.w_max];
    zvfull = [cba_constraints.zv; zeros(nr,1); zeros(nm,1)];

    
    %% inequality cba_constraints: 
    %% positive objective
    %% external signs
    
    A = [zv', zeros(1,nr+nm) ; ...
         double([-diag(cba_constraints.ext_sign(ind_extsigns)) * N(ind_extsigns,:), ...
                 zeros(length(ind_extsigns),nr+nm)])];
    
    b = [0; zeros(length(ind_extsigns),1)];
    
    %% inequality cba_constraints: 
    %% positive objective
    %% external signs
    
    A = [-zv', zeros(1,nr+nm) ; ...
         double([-diag(cba_constraints.ext_sign(ind_extsigns)) * N(ind_extsigns,:), ...
                 zeros(length(ind_extsigns),nr+nm)])];
    
    b = [0; zeros(length(ind_extsigns),1)];
    
end


% ----------------------------------------------------------
% optimisation

ppp.Ntrans  = full(N');
ppp.ind_v   = ind_v;
ppp.ind_dmu = ind_dmu;
ppp.ind_y   = ind_y;
ppp.ind_w   = ind_w;
ppp.zv      = zv;

opt = optimset('MaxFunEvals',10000);

switch cba_options.objective
  
  case 'fba',

    switch cba_options.cba_conditions
      case 'y', %% use constraint function cba2f
        x0 = rand(3*nr,1)-.5;
        weights = regularisation_weight * ones(size(x0));
        ppp.epsilon = -0.1;
        [xtry,negvalue,exitflag] = fmincon(@(xx)cbaf1(xx,zvfull,weights), x0, A, b, Aeq, beq, xmin, xmax, @(xx)cbaf2(xx,ppp), opt);
        %% run, with slightly relaxed feasibility criterion
        ppp.epsilon = -0.001;
        [x,negvalue,exitflag] = fmincon(@(x)cbaf1(x,zvfull,weights), xtry, A, b, Aeq, beq, xmin, xmax, @(xx)cbaf2(xx,ppp), opt);
      case 'w', %% use constraint function cba2fw
        x0 = rand(2*nr+nm,1)-.5;
        weights = regularisation_weight * ones(size(x0));
        ppp.epsilon = -0.1;
        [xtry,negvalue,exitflag] = fmincon(@(xx)cbaf1(xx,zvfull,weights), x0, A, b, Aeq, beq, xmin, xmax, @(xx)cbaf2w(xx,ppp), opt);
        %% run, with slightly relaxed feasibility criterion
        ppp.epsilon = -0.001;
        [x,negvalue,exitflag] = fmincon(@(x)cbaf1(x,zvfull,weights), xtry, A, b, Aeq, beq, xmin, xmax, @(xx)cbaf2w(xx,ppp), opt);
    end
    
  case 'fit',

    v_mean = cba_constraints.v_mean;
    v_std  = cba_constraints.v_std;
    switch cba_options.cba_conditions
      case 'y', %% use constraint function cba2f
        x0 = rand(3*nr,1)-.5;
        ppp.epsilon = -0.1;
        [xtry,value,exitflag] = fmincon(@(xx)cbaf3(xx,v_mean,v_std), x0, A, b, Aeq, beq, xmin, xmax, @(xx)cbaf2(xx,ppp), opt);
        %% run, with slightly relaxed feasibility criterion
        ppp.epsilon = -0.001;
        [x,value,exitflag] = fmincon(@(x)cbaf3(x,v_mean,v_std), xtry, A, b, Aeq, beq, xmin, xmax, @(xx)cbaf2(xx,ppp), opt);
      case 'w', %% use constraint function cba2fw
        x0 = rand(2*nr+nm,1)-.5;
        weights = regularisation_weight * ones(size(x0));
        ppp.epsilon = -0.1;
        [xtry,value,exitflag] = fmincon(@(xx)cbaf3(xx,v_mean,v_std), x0, A, b, Aeq, beq, xmin, xmax, @(xx)cbaf2w(xx,ppp), opt)
        %% run, with slightly relaxed feasibility criterion
        ppp.epsilon = -0.001;
        [x,value,exitflag] = fmincon(@(x)cbaf3(x,v_mean,v_std), xtry, A, b, Aeq, beq, xmin, xmax, @(xx)cbaf2w(xx,ppp), opt)
    end
    
end

% ----------------------------------------------------------
% format the output

v     = nan * ones(nr,1);
dmu   = nan * ones(nr,1);
y     = nan * ones(nr,1);
w     = nan * ones(nm,1);
value = nan;

if exitflag<=0,
  display('CBA failed');
else,
  v   = x(ind_v);
  dmu = x(ind_dmu);
  switch cba_options.objective
    case 'fba',
      value = -negvalue; 
  end
  switch cba_options.cba_conditions,
    case 'y', y = x(ind_y);
    case 'w', w = x(ind_w);
  end

  if value < 0, warning('Objective function is negative'); end
  if sum(v.*dmu > 0), warning('Solution violates EBA constraint'); end
  switch cba_options.cba_conditions,
    case 'y', 
      if sum(v.*y < 0), warning('Solution violates CBA constraint'); end
    case 'w',
      if sum(v.*[N'*w+zv] < 0), warning('Solution violates CBA constraint'); end
  end
  
end

result.value = value; 
result.v     = v;
result.dmu   = dmu;

switch cba_options.cba_conditions,
  case 'y',    result.y  = y;
  case 'w',    result.w  = w;
end
