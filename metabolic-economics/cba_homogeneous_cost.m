function [w, delta_w, y, zx] = cba_homogeneous_cost(network,v,cba_constraints,y_given,method)

% CBA_HOMOGENEOUS_COST - Determine economic potentials from the principle of even cost
%
% [w, delta_w, y, zx] = cba_homogeneous_cost(network, v, cba_constraints, y_given, method)
%
% Assumptions: 
%  o Balanced distribution of enzyme cost (min = sum(u^2))
%  o Benefit arises only from fluxes and external production (not from concentrations)
%
% Input
%   y_given: (optional) scalar or vector (#reactions x 1) of enzyme costs to be approximated
%
%   method: {'cost', 'force'}  if method == 'force', the vector y_given does not represent
%                                    the cost y, but the force y./v; in this case, it has to 
%                                    be explicitly given as a vector

% THE LATTER OPTION HAS NOT BEEN TESTED YET


eval(default('y_given','[]','method','''cost'''));

ind_act = find(v ~=0 );
ind_int = find(network.external ==0);
ind_ext = find(network.external ==1);

cba_constraints = cba_update_constraints(cba_constraints,network.N(find(network.external),:),network);

switch method, 
  case 'cost',
    %% if cost unknown -> choose according to benefit function
    if isempty(y_given), y_given = [cba_constraints.zv' * v ] / length(ind_act); end

    %% cost: scalar ->  vector 
    if length(y_given) == 1, y_given = y_given * double(v~=0); end

  case 'force',
    y_given(ind_act) = y_given(ind_act) ./ v(ind_act);
end

if sum(y_given(ind_act) <= 0), error('Impossible enzyme costs'); end 

% rescale benefit function to match the predefined costs

cba_constraints.z_ext = cba_constraints.z_ext * [sum(y_given)] / [cba_constraints.zv'*v];
cba_constraints.z_int = cba_constraints.z_int * [sum(y_given)] / [cba_constraints.zv'*v];
cba_constraints.zv    = cba_constraints.zv    * [sum(y_given)] / [cba_constraints.zv'*v];

v_act        = v(ind_act);
y_given_act  = y_given(ind_act);
N_int        = network.N(ind_int,ind_act);


% quadratic function min = wc' * A * wc + a' * wc
% s.t.                 B * wc <= b;
%(without constraints the solution would be simply ww_int = - pinv(A) * a;)

epsilon = 10^-3;

A  = N_int * diag( [v_act.^2] ./ y_given_act ) * N_int';
a  = N_int * diag( [v_act.^2] ./ y_given_act ) * cba_constraints.zv(ind_act);
B  = - [diag(v_act) * N_int'];
nr = length(v_act);
b  = - [epsilon * ones(nr,1) - diag(v_act) * cba_constraints.zv(ind_act)];

% regularisation

if find(eig(A)==0),
  alpha = 10^-8; A = A + alpha * eye(size(A));
end

%  if exist('cplexqp','file'),
%    opt = cplexoptimset('Display','off','Algorithm','interior-point-convex');
%    w_int = cplexqp(A,a,[],[],[],[],[],[],[],opt); 
%  else,
    opt = optimset('Display','off','Algorithm','interior-point-convex');
    w_int = quadprog(A,a,[],[],[],[],[],[],[],opt); 
%  end

if(sum(double(B * w_int > b))),
  %% constraints violated
  if exist('cplexqp','file'),
    opt = cplexoptimset('Algorithm', 'interior-point-convex','Display','off');
    w_int = cplexqp(A,a,B,b,[],[],[],[],[],opt);
  else
    opt = optimset('Algorithm', 'interior-point-convex','Display','off');
    w_int = quadprog(A,a,B,b,[],[],[],[],[],opt);
  end
end


% ---------------------

w(ind_ext,1) = cba_constraints.z_ext;
w(ind_int)   = w_int;

delta_w = network.N'*w;

y           = v .* [delta_w + cba_constraints.z_int];
zx          = zeros(size(w)); 
zx(ind_ext) = cba_constraints.z_ext;

% ---------------------

if find(y(v~=0) == 0), warning('Zero enzyme cost encountered'); end
