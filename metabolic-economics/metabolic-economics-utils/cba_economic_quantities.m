function res = cba_economic_quantities(network, c_opt, u_opt, cba_constraints)

% CBA_METABOLITE_AND_ENZYME_PRICES - Convenience function for computing economic quantities
%
% res = cba_economic_quantities(network,c_opt,u_opt,cba_constraints)
%
% convenience function, calls cba_economic_potentials and the respective cost function

u_opt(find(u_opt/max(u_opt)<10^-5)) = 0;

res.u_opt          = u_opt;
network.kinetics.u = u_opt;

[nm,nr] = size(network.N);

if isempty(cba_constraints.zc),
  cba_constraints.zc = zeros(nm,1);
end

[res.wc, v_opt, res.dw_du, dv_du] = cba_economic_potentials(network, c_opt, cba_constraints.zv, cba_constraints.zc, cba_constraints.ind_controllable);

res.v_opt  = v_opt;
res.f_bene = cba_constraints.zv' * v_opt;

[res.f_cost, res.f_cost_gradient, res.f_cost_hessian] = feval(cba_constraints.fitness_options.cost_function,u_opt,cba_constraints.fitness_options);

res.h_hat = res.f_cost_gradient .* u_opt;
res.y     = res.f_cost_gradient .* u_opt ./ v_opt;

res.w = res.wc;
res.w(find(network.external)) = res.w(find(network.external)) + cba_constraints.z_ext;