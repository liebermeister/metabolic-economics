function [res,res_eq] = cba_opt_constraints(u,c_init,network,constraints)

global my_c

network.kinetics.u = u;
[s,v,sdot] = network_steady_state(network,my_c);

res_sign   = - v(isfinite(constraints.v_sign)) .* constraints.v_sign(isfinite(constraints.v_sign));
res_min    =  constraints.v_min(isfinite(constraints.v_min)) - v(isfinite(constraints.v_min));

res = [res_sign; res_min];

res_eq = [];