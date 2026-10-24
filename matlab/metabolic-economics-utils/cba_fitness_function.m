function [f, s, v, fplus, fminus, zx] = cba_fitness_function(u, c_init, network, constraints)

% [f,s,v,fplus,fminus,zx, s, v] = cba_fitness_function(u,c_init,network,constraints)
%
% Fitness function used in CBA. Parameters are read from the fields of constraints

global my_c
global my_v

c_init = my_c;

network.kinetics.u = u;
[s, v, sdot]       = network_steady_state(network, c_init);
fplus              = constraints.zv' * v;
if isfield(constraints,'fitness_options'),
  fminus = feval(constraints.fitness_options.cost_function,u,constraints.fitness_options);
else
  fminus = constraints.fx' * sum(u) + 0.5 * constraints.fxx * sum(u.^2);
end
f                  = fplus - fminus;

% provide gradient?
if 0,
  R   = basic_control_analysis(network,s);
  RJu = R.RJ(:,1:length(u));
  f_plus_grad  = [constraints.zv' * RJu]';
  f_minus_grad = constraints.fx * ones(size(u)) + constraints.fxx * u;
  f_grad       = f_plus_grad - f_minus_grad;
end

if nargout > 5,
  nr = length(u);
  zx = constraints.fx' * ones(nr,1) + constraints.fxx * u;
end

% punish solutions violating stationarity
epsilon = 10^-3;
error   = max(abs(sdot ./ [s(find(1-network.external)) + epsilon]));

% v
% sdot
% s(find(1-network.external))

if error > 0.0001, 
  f = -10^8*error; 
else
  my_c = s;
  my_v = v;
end

display(sprintf('f: %f f_bene: %f f_cost: %f',f,fplus, fminus))
if rand < 0.2,
  gp = struct('arrowvalues',v,'actstyle','fixed','flag_edges',1);
  netgraph_concentrations(network,s,u,1,gp); drawnow
end
