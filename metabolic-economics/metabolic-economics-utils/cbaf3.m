function f = cbaf3(x,v_mean,v_std)

% constraints for EBA and CBA
% f = cbaf3(x,v_mean,v_std)

ind = find(isfinite(v_mean));

f = sum([[x(ind)-v_mean(ind)]./v_std(ind)].^2);
