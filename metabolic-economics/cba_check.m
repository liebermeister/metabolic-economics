function isFeasible = cba_check(N, external, cba_constraints, cba_options, v, dmu, hv)

% CBA_CHECK - Check flux mode, chem. pot. differences, and spec. flux costs for feasibility
%
% isFeasible = cba_check(N, external, cba_constraints, cba_options, v, dmu, hv)
%
% Check whether flux mode v, chemical potentials dmu, and specific flux costs hv
% satisfy the thermodynamic and economic constraints 

% !!!
% care about zero fluxes ???

%  ------------------------------------------------------------
% initialisation

%rand('state',cba_options.seed);

[nm,nr]      = size(N);

if ~isfield(cba_constraints,'dmu_min'), cba_constraints.dmu_min = -10 *ones(nr,1); end
if ~isfield(cba_constraints,'dmu_max'), cba_constraints.dmu_max =  10 *ones(nr,1);; end
if ~isfield(cba_constraints,'y_min'),   cba_constraints.y_min   = -10 *ones(nr,1);; end
if ~isfield(cba_constraints,'y_max'),   cba_constraints.y_max   =  10 *ones(nr,1);; end

cba_constraints.v_min(find(cba_constraints.v_sign== 1)) = 0;
cba_constraints.v_max(find(cba_constraints.v_sign==-1)) = 0;

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

weights      = ones(3*nr,1);     % weights in regularisation term


% ------------------------------------------------------------
% vector x contains all v, dmu, and hv
% indices within vector:

ind_v   =          1:nr;
ind_dmu =    nr + (1:nr);
ind_hv  = 2* nr + (1:nr);

% ------------------------------------------------------
% cba_constraints:
%
% limits           : cba_constraints.v_min   <=  x(ind_v)   <= cba_constraints.v_max
%                  : cba_constraints.dmu_min <=  x(ind_dmu) <= cba_constraints.dmu_max
%                  : cba_constraints.y_min   <=  x(ind_hv)  <= cba_constraints.y_max
%
% fixed values     : x(ind_v_fix)         = cba_constraints.v_fix(ind_v_fix)
% stationary fluxes: N_int * x(ind_v )    = 0
% dmu constraint   : Ktot' * x(ind_hv)    = 0
% hv constraint    : Kint' * x(ind_hv)    = Kint' * zv
%
% feasible signs   : diag(x(ind_v)) * N' * x(ind_dmu) <= - epsilon < 0   where v == 0
%                  : diag(x(ind_v)) * N' * x(ind_hv)  >=   epsilon > 0   where v == 0
%
% regularisation   : min != sum(x.^2);


% ------------------------------------------------------
% rewrite this in the form 
% 
% Aeq * x  = beq
%   A * x <= b
%       C  = cbaf2(x,N)   (for signs)
%     min != cbaf1(x)     (regularisation)


xmin = [cba_constraints.v_min; cba_constraints.dmu_min; cba_constraints.y_min];
xmax = [cba_constraints.v_max; cba_constraints.dmu_max; cba_constraints.y_max];

Aeq  = double([my_eye(ind_v_fix,:),  zeros(length(ind_v_fix),2*nr)      ; ...
               N_int,                 zeros(n_int,2*nr)                ; ...
               zeros(nKtot,nr),       Ktot',          zeros(nKtot,nr)  ; ...
               zeros(nKint,2*nr),                     Kint'             ]);

beq  = double([cba_constraints.v_fix(ind_v_fix);... 
               zeros(n_int,1)           ;...
               zeros(nKtot,1)    ;...
               Kint' * zv ]);

A = double([-diag(cba_constraints.ext_sign(ind_extsigns)) * N(ind_extsigns,:), ...
            zeros(length(ind_extsigns),2*nr)]);

b = zeros(length(ind_extsigns),1);


% ------------------------------------------------------------------------

x = [v; dmu; hv];

isFeasible = 1;

ppp.N = full(N);
ppp.ind_v   = ind_v  ;
ppp.ind_dmu = ind_dmu;
ppp.ind_hv  = ind_hv ;
ppp.epsilon = 0;

if sum(x < xmin) + sum(x > xmax),          warning('Limit constraints violated');      isFeasible = 0; end 
if norm( Aeq * x - beq)/length(beq)>10^-5, warning('Equality constraints violated');   isFeasible = 0; end 
if sum(A * x > b),                         warning('Inequality constraints violated'); isFeasible = 0; end 
if sum(cbaf2(x,ppp) > 0),                  warning('Sign constraints violated');       isFeasible = 0; end 
