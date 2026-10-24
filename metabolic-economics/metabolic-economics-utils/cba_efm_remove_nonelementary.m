function M = efm_remove_nonelementary(M,E)

% M = efm_remove_nonelementary(M,E)

if length(M) * length(E),
  nonelementary = zeros(1,size(M,2));
  for it = 1:size(E,2),
    nonelementary = nonelementary + prod(M(find(E(:,it)),:),1);
  end
  M = M(:,nonelementary==0);
end
