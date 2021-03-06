function [C, Ceq] = cbaf2(x,p)

% constraints for EBA and CBA 
% [C, Ceq] = cbaf2(x,p)

C = [ x(p.ind_v) .* x(p.ind_dmu) + p.epsilon; ...
    - x(p.ind_v) .* x(p.ind_y)   + p.epsilon ];

% no EBA constraints where v = 0

C(find(x(p.ind_v)==0)) = -1;

% no CBA constraints where v = 0

C(find(x(p.ind_v)==0) + length(p.ind_v))=-1;

Ceq     = [];
