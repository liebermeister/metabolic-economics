function cba_constraints = cba_update_constraints(cba_constraints,Next,network)

% cba_constraints = cba_update_constraints(cba_constraints,Next)

cba_constraints.dmu_min(cba_constraints.dmu_sign>0) = 0;
cba_constraints.dmu_max(cba_constraints.dmu_sign<0) = 0;

if length(cba_constraints.z_ext),
  zv = Next' * cba_constraints.z_ext + cba_constraints.z_int;
else,
  cba_constraints.z_ext = zeros(size(Next,1),1);
  zv = cba_constraints.z_int;
end  

if isfield(cba_constraints,'zv'),
  if length(cba_constraints.zv),
    if sum(zv ~= cba_constraints.zv)
      warning('Changing existing entry zv');
    end
  end
end

cba_constraints.zv = zv;

cba_constraints = fba_update_constraints(cba_constraints,network);

