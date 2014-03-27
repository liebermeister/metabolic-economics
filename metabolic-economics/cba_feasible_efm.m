function [isFeasible, ind_conflict] = cba_feasible_efm(v, C_efm)

% CBA_FEASIBLE_EFM - Test flux mode for economic feasibility (given non-beneficial flux modes)
%
% [isFeasible, ind_conflict] = cba_feasible_efm(v, C_efm)
%
% Test a flux mode v for economic feasibility by comparing it to
% elementary non-beneficial flux modes (in matrix C_efm)
%
% ind_conflict: indices of those modes in C_efm that are in conflict with v

if find(v==0), 

  warning('Inactive reactions encountered. Feasibility cannot be directly tested');
  isFeasible   = nan; 
  ind_conflict = [];

else,

  isFeasible   = 1;
  ind_conflict = [];

  for it = size(C_efm,2),
    ind_intersect = find( [v~=0] .* [C_efm(:,it)~=0]);
    if length(unique(sign(v(ind_intersect).*C_efm(ind_intersect,it))))==1,
      isFeasible   = 0;
      ind_conflict = [ind_conflict, it];
    end
  end
  
end
