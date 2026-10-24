function M = efm_remove_duplicates(M)

% EFM_REMOVE_DUPLICATES - Remove duplicate columns in Elementary Flux Mode Matrix
%
% M = efm_remove_duplicates(M)

if length(M),
  
  duplicate = zeros(size(M,2),1);

  %% trick!!
  zeta = randn(1,size(M,1)) * M;
  [dum,order] = sort(zeta);
  M = M(:,order);
  
  n_check = 20; % how many duplicates of one vector are expected at most?

  for it = 2:size(M,2),
    ind_check = max(1,it-100):it-1;
    %%  duplicate(it) = M(:,ind_check)'*M(:,it) == M(:,it)'*M(:,it);
    duplicate(it) = sum( full( sum( abs( M(:,ind_check) - repmat(M(:,it),1,length(ind_check))) ,1))==0)>0;
  end
  M = M(:,duplicate==0);
  
end
