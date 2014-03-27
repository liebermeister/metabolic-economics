 function [g, geq, Ec_un, beta_M, beta_A, beta_I, c, mu, CJ] = cba_flux2model_mrl_constraint_xi(xi,Mplus,Mminus,Wplus,Wminus,v,N,ind_ext,cba_constraints,kinetic_law,h)

% function [g, geq, Ec_un, beta_M, beta_A, beta_I, c, mu, CJ] = cba_flux2model_mrl_constraint_xi(xi,Mplus,Mminus,Wplus,Wminus,v,N,ind_ext,cba_constraints,kinetic_law)
%
% choose state parameters (mu0, c, alphas) such that
%  o all flux control coefficient on benefit are positive
%  o all reaction affinities and fluxes have the same signs

eval(default('h','ones(size(v))'));

[nr, nm] = size(Mplus);

M       = Mplus + Mminus;
alpha_M = zeros(size(M)); 
alpha_A = zeros(size(M)); 
alpha_I = zeros(size(M)); 
beta_M  = zeros(size(M)); 
beta_A  = zeros(size(M)); 
beta_I  = zeros(size(M)); 

ind_KM = find(M);
ind_KA = find(Wplus);
ind_KI = find(Wminus);
n_KM   = length(ind_KM);
n_KA   = length(ind_KA);
n_KI   = length(ind_KI);

mu0  = xi(1:nm);
ln_c = xi(nm+[1:nm]);
c    = exp(ln_c);

mu   = mu0 + RT*ln_c;
A    = -N'*mu;
A(A==0) = 10^-10; % simple fix to avoid divergence
zeta = exp(h.*A/RT);

E_sc_T = diag(zeta./[zeta-1]) * Mplus - diag(1./[zeta-1]) * Mminus;

alpha_M(ind_KM) = 1-xi(2*nm+[1:n_KM]);
beta_M( ind_KM) =   xi(2*nm+[1:n_KM]);

alpha_A(ind_KA) = 1-xi(2*nm+n_KM + [1:n_KA]);
beta_A( ind_KA) =   xi(2*nm+n_KM + [1:n_KA]);

alpha_I(ind_KI) = 1-xi(2*nm+n_KM + n_KA + [1:n_KI]);
beta_I( ind_KI) =   xi(2*nm+n_KM + n_KA + [1:n_KI]);

% -----------------------------------------------------------
% copy from compute_modular_elasticities

switch kinetic_law,  

  case {'cs','ms'}
    alpha_M( find((Mplus==0).*(Mminus==0)) ) = 1;  % irrelevant elements
    psi_plus    = prod((1./alpha_M) .^ Mplus, 2);
    psi_minus   = prod((1./alpha_M) .^ Mminus,2);

  case {'ds','fd'},
    alpha_M( find((Mplus==0).*(Mminus==0)) ) = 1/2; % irrelevant elements
    xi_plus  = prod((1./alpha_M-1) .^ Mplus, 2);
    xi_minus = prod((1./alpha_M-1) .^ Mminus,2);

  otherwise error('Unknown kinetic law');

end

switch kinetic_law,  
  case 'cs',
    D               = psi_plus + psi_minus - 1;
    E_sc_Dc_kinetic = diag(1./D) * beta_M .* [diag(psi_plus) * Mplus + diag(psi_minus) * Mminus];
  case 'ms',
    D               = psi_plus .* psi_minus;
    E_sc_Dc_kinetic = beta_M   .* [ Mplus + Mminus ];
  case 'ds',
    D               = xi_plus + xi_minus + 1;
    E_sc_Dc_kinetic = diag(1./D) * [ diag(xi_plus) * Mplus + diag(xi_minus) * Mminus ];
  case 'rp',
    D               = ones(nm,1);
    E_sc_Dc_kinetic = zeros(nr,nm);
  case 'fd',
    D               = sqrt(xi_plus .* xi_minus);
    E_sc_Dc_kinetic = 1/2 * [ Mplus + Mminus ];
end

E_sc_c_regulation = alpha_A .* Wplus - beta_I .* Wminus;

% -----------------------------------------------------------

Ec_un = diag(v) * [E_sc_T - E_sc_Dc_kinetic + E_sc_c_regulation] * diag(1./c);

ind_int = setdiff(1:nr,ind_ext);
N_int   = N(ind_int,:);
E_int   = Ec_un(:,ind_int);

CJ      = eye(nr) - E_int * pinv(N_int*E_int) * N_int;

g       = [- diag(sign(v)) * CJ' * cba_constraints.zv + cba_constraints.zx_scaled_min; ...
           ];%-sign(v .* A) .* double(v~=0) - 0.1 ];

geq     = [];
