function [g, geq] = cba_flux2model_constraint(E_vector,N,ind_ext,zv,epsilon)

% [g, geq] = cba_flux2model_constraint(E_vector,N,ind_ext,zv,epsilon)

nr   = length(zv);
E    = -N';
E(find(E)) = E_vector;

ind_int = setdiff(1:nr,ind_ext);
N_int = N(ind_int,:);
E_int = E(:,ind_int);
CJ    = eye(nr) - E_int * inv(N_int*E_int) * N_int;

g     = - CJ' * zv + epsilon;
geq   = [];
