function [w, v, dw_du,  dv_du, c] = cba_economic_potentials(network, c, zv, zc, ind_controllable)

% CBA_ECONOMIC_POTENTIALS - Compute the economic potentials and their sensitivities in a kinetic model
%
% [w, v, dw_du,  dv_du] = cba_economic_potentials(network, c, cba_constraints,ind_controllable)
%
% w contains only the internal economic potentials
% Compute the economic potentials for a kinetic model in a given stationary state
% as well as their derivatives dw_du with respect to the controllable enzyme levels
% 
% The calculation works only for linear benefit functions (defined by weights zv, zc)


% ---------------------------------------------------------
% to be sure ... put system into steady state

[c,v]   = network_steady_state(network,c,200);
[nm,nr] = size(network.N);
ind_int = find(network.external==0);


% ---------------------------------------------------------
% Compute economic potentials w
% phi: virtual supply fluxes for independent metabolites 


Ec = elasticities(network,c);
[CJ, CS, L_int, NR_int, M, Madj, indp_among_internal] = control_coefficients(network.N, Ec, network.external);

CS_int_phi = - L_int * pinv(M);
CS_phi     = zeros(nm,length(indp_among_internal));
CS_phi(ind_int,:) = CS_int_phi;
CJ_phi     = Ec(:,ind_int) * CS_int_phi;
w_ind      = [zv' * CJ_phi + zc' * CS_phi]';
w          = zeros(nm,1);
w(ind_int(indp_among_internal)) = w_ind;


% ---------------------------------------------------------
% Compute sensitivities dw_du and dv_du

if nargout > 2,

  %% NUMERICAL DERIVATIVES
  uu = network.kinetics.u;
  dd = 10^-8;

  for itt = 1:nr,
    my_du1        = 0 * uu;
    my_du1(itt)   =  uu(itt) * dd;
    my_du2        = -my_du1;
    my_du1(uu==0) =  dd;
    my_du2(uu==0) =   0;
    network.kinetics.u = uu + my_du1;
    [my_w1, my_v1] = cba_economic_potentials(network, c, zv, zc, ind_controllable);
    network.kinetics.u = uu + my_du2;
    [my_w2, my_v2] = cba_economic_potentials(network, c, zv, zc, ind_controllable);
    dw_du(:,itt) = [my_w1 - my_w2] / [my_du1(itt) - my_du2(itt)];
    dv_du(:,itt) = [my_v1 - my_v2] / [my_du1(itt) - my_du2(itt)];
  end

  %% ANALYTICAL CALCULATION:
  %% R is not actually control matrix, but "response" to enzyme levels
  % Ec                         = elasticities(network,c);
  % [Ec,Ep,parameters,Ecc,Ecp,Epp,p] = elasticities(network,c);
  % Ecu = Ecp(:,:,1:nr);
  % v = network_velocities(c,network);
  % rho_2 = diag(1./v) * Ec;
  %% VORSICHT: GAMMA-FORMEL IST eventuell noch  FALSCH!!!!
  %% es wird keine term berücksichtigt, der den direkten einfluss
  %von phi mit dem indirekten von u koppelt. nochmal nachrechnen!
  % Gamma =   tensor_product(tensor_product(Ecc,CS_phi),CS,2,1) ...
  %         + permute(tensor_product(Ecu,CS_phi,2,1),[1,3,2]); 
  % RS_phi_u = tensor_product(CS,Gamma);
  % RJ_phi_u = tensor_product(Ec,RS_phi_u);
  % 
  % dw_du_ind = squeeze(tensor_product(zv', RJ_phi_u) + tensor_product(zc', RS_phi_u) );
  % dw_du = zeros(nm,nr);
  % dw_du(ind_int(indp_among_internal),:) = dw_du_ind;

end
