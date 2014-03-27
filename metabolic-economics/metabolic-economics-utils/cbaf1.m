function f = cbaf1(x,zv,weights)

% constraints for EBA and CBA
% f = cbaf1(x,zv,weights)

f = zv'*x + sum(weights .* x .^ 2);
