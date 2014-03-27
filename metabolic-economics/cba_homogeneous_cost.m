function [w, delta_w, q, zx] = cba_homogeneous_cost(network,v,cba_constraints,q_given,method)

% CBA_HOMOGENEOUS_COST - Determine economic potentials from the principle of even cost
%
% [w, delta_w, q, zx] = cba_homogeneous_cost(network, v, cba_constraints, q_given, method)
%
% Assumptions: 
%  o Balanced distribution of enzyme cost (min = sum(u^2))
%  o Benefit arises only from fluxes and external production (not from concentrations)
%
% Input
%   q_given: (optional) scalar or vector (#reactions x 1) of enzyme costs to be approximated
%
%   method: {'cost', 'force'}  if method == 'force', the vector q_given does not represent
%                                    the cost q, but the force q./v; in this case, it has to 
%                                    be explicitly given as a vector

% THE LATTER OPTION HAS NOT BEEN TESTED YET


eval(default('q_given','[]','method','''cost'''));

ind_act = find(v ~=0 );
ind_int = find(network.external ==0);
ind_ext = find(network.external ==1);

cba_constraints = cba_update_constraints(cba_constraints,network.N(find(network.external),:));

switch method, 
  case 'cost',
    %% if cost unknown -> choose according to benefit function
    if isempty(q_given), q_given = [cba_constraints.zv' * v ] / length(ind_act); end

    %% cost: scalar ->  vector 
    if length(q_given) == 1, q_given = q_given * double(v~=0); end

  case 'force',
    q_given(ind_act) = q_given(ind_act) ./ v(ind_act);
end

if sum(q_given(ind_act) <= 0), error('Impossible enzyme costs'); end 

% rescale benefit function to match the predefined costs

cba_constraints.z_ext = cba_constraints.z_ext * [sum(q_given)] / [cba_constraints.zv'*v];
cba_constraints.z_int = cba_constraints.z_int * [sum(q_given)] / [cba_constraints.zv'*v];
cba_constraints.zv    = cba_constraints.zv    * [sum(q_given)] / [cba_constraints.zv'*v];

v_act        = v(ind_act);
q_given_act  = q_given(ind_act);
N_int        = network.N(ind_int,ind_act);


% quadratic function min = wc' * A * wc + a' * wc
% s.t.                 B * wc <= b;
%(without constraints the solution would be simply ww_int = - pinv(A) * a;)

A = N_int * diag( [v_act.^2] ./ q_given_act ) * N_int';
a = N_int * diag( [v_act.^2] ./ q_given_act ) * cba_constraints.zv(ind_act);
B = - [diag(v_act) * N_int'];
epsilon = 10^-3;
nr = length(v_act);
b = - [epsilon * ones(nr,1) - diag(v_act) * cba_constraints.zv(ind_act)];

% regularisation

if find(eig(A)==0),
  alpha = 10^-8; A = A + alpha * eye(size(A));
end

opt = optimset('Display','off','Algorithm','interior-point-convex');
w_int = quadprog(A,a,[],[],[],[],[],[],[],opt); 

if(sum(double(B * w_int > b))),
  %% constraints violataed
  opt = optimset('Algorithm', 'interior-point-convex','Display','off');
  w_int = quadprog(A,a,B,b,[],[],[],[],[],opt);
end


% ---------------------

w(ind_ext,1) = cba_constraints.z_ext;
w(ind_int)   = w_int;

delta_w = network.N'*w;

q           = v .* [delta_w + cba_constraints.z_int];
zx          = zeros(size(w)); 
zx(ind_ext) = cba_constraints.z_ext;

% ---------------------

if find(q(v~=0) == 0), warning('Zero enzyme cost encountered'); end
 
