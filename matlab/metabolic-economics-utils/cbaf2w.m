function [C, Ceq] = cbaf2w(x,p)

% constraints for EBA and CBA
% [C, Ceq] = cbaf2w(x,p)

C = [ x(p.ind_v) .* x(p.ind_dmu) + p.epsilon; ...
    - x(p.ind_v) .* [p.Ntrans * x(p.ind_w) + p.zv] + p.epsilon ];

% no EBA constraints where v = 0

C(find(x(p.ind_v)==0)) = -mean(abs(C));

% no CBA constraints where v = 0

C(find(x(p.ind_v)==0) + length(p.ind_v))=-mean(abs(C));

Ceq = [];
