function [isFeasible, res] = cba_check_model(network, c, u, v, cba_constraints)

% CBA_CHECK_MODEL - Check kinetic model for being in an enzyme-optimal state
%
% [isFeasible,res] = cba_check_model(network,c,u,v,cba_constraints)
%
% Check a kinetic model for a feasible stationary state
%
% 'res' contains various kinds of information about the state: 
%
%   res.c_adjusted
%   res.v_adjusted
%   res.u_adjusted
%   res.c_mismatch_norm
%   res.v_mismatch_norm
%   res.mismatch_v_constraints
%   res.stationarity_mismatch_norm
%   res.flux_sign_mistakes
%   res.dynamically_instable
%   res.negative_control
%   res.R_f_u
%   res.economically_unstable
%   res.fitness_gradient
%   res.fitness_hessian
%   res.benefit_hessian
%   res.cost_hessian


% stationarity

% set small enzyme levels to zero

u(u/median(u) < 0.001) = 0;

[nm,nr] = size(network.N);

network.kinetics.u = u;

[c_st,v_st]    = network_steady_state(network,c);

res.c_adjusted = c_st;
res.v_adjusted = v_st;
res.u_adjusted = u;

res.c_mismatch_norm = norm(c_st - c)/norm(c);
res.v_mismatch_norm = norm(v_st - v)/norm(v);

res.mismatch_v_constraints     = sum(v<cba_constraints.v_min) + sum(v>cba_constraints.v_max);
res.stationarity_mismatch_norm = norm(network.N(network.external==0,:) * v_st);
res.flux_sign_mistakes         = find([v_st~=0] .*  [sign(v_st .* [ log(network.kinetics.Keq) - network.N'*log(c_st)]) ==-1]);

R   = basic_control_analysis(network,c_st);
RJu = R.RJ(:,1:nr);

res.dynamically_instable =  sum(real(eig(R.M))>0)>0;

R_f_u =  cba_constraints.zv' * RJu;
R_f_u(abs(R_f_u)<0.001*median(abs(R_f_u))) = 0;

res.negative_control = find([u~=0].*R_f_u' <0);

if res.negative_control,
  res.R_f_u = R_f_u;
end

if length(res.negative_control),
  figure(10); netgraph_concentrations(network,[], [u~=0].*R_f_u',1); title('Fitness/enzyme response coefficients');
end

[f_cost, grad_cost, H_cost] = feval(cba_constraints.fitness_options.cost_function, u, cba_constraints.fitness_options);

fitness_gradient = [cba_constraints.zv' * RJu] - grad_cost';
res.fitness_gradient_active_norm = norm(fitness_gradient .* [u ~=0]');

[Ec,Ep,parameters,Ecc,Ecp,Epp,p] = elasticities(network,c);
[CJ, CS] = control_coefficients(network.N, Ec, network.external);
[RS,RJ,RS2,RJ2] = response_coefficients(CS,Ec,Ep,Ecc,Ecp,Epp);

H_benefit  = squeeze(tensor_product(cba_constraints.zv', RJ2(:,1:nr,1:nr)));
H_fitness  = H_benefit - H_cost;

res.economically_unstable = sum(eig(H_fitness(find(u>0),find(u>0)))>0)>0;

res.fitness_gradient = fitness_gradient;
res.fitness_hessian  = H_fitness;
res.benefit_hessian  = H_benefit;
res.cost_hessian     = H_cost;

isFeasible = ...
    [res.c_mismatch_norm^2/nm <10^-5] ...
    * [res.v_mismatch_norm^2/nr <10^-5] ...
    * [res.stationarity_mismatch_norm^2/nr <10^-5] ...
    * [res.mismatch_v_constraints == 0] ...
    * [length(res.flux_sign_mistakes) == 0] ...
    * [res.dynamically_instable == 0] ...
    * [length(res.negative_control) == 0] ...
    * [res.fitness_gradient_active_norm^2/sum(u>0) < 10^-5] ...
    * [res.economically_unstable == 0];
