% Compute the clustering error.
%
%-------------------------------------------------------
% Input:
%       A:  predicted cluster assignments
%       A0 : N true cluster assignments
%       truth:  truth cluster indicators
%     
%
% Output:
%       CE:  clustering error
%
function CE = compute_CE(A,A0)

if size(A,2) == 1
    A = A';
end

if size(A0,2) == 1
    A0 = A0';
end

L  = max(A);
L0 = max(A0);

N = numel(A);

Aerror = zeros(L,max(L,L0));
for i=1:L
    for j=1:L0
        Aerror(i,j) = nnz(A0(A == i) ~= j);
    end
    Aerror(i,L0+1:end) = nnz(A == i);
end

if max(L,L0) <= 10
    % for <=10 labels, compute error for every permutation
    perm = perms(1:max(L,L0));
    perm = perm(:,1:L);
    
    ind_set = repmat((0:L-1)*L,size(perm,1),1) + perm;
    [CE,~]  = min(sum(Aerror(ind_set),2));  % Find the best permutation of label indices that minimize the disagreement
    CE = CE / N;
    
else
    % for >10 labels, compute approximate error using hill climbing
    % assumed L=L0
    swap = [];
    for i=2:L
        swap = [swap; repmat(i,i-1,1) (1:i-1)'];
    end
    
    CE = N;
    for i=1:10
        perm = randperm(L);
        
        for j=1:1e5
            [m,ind] = min(Aerror(sub2ind(size(Aerror),swap(:,1),perm(swap(:,2))')) ...
                + Aerror(sub2ind(size(Aerror),swap(:,2),perm(swap(:,1))')) ...
                - Aerror(sub2ind(size(Aerror),swap(:,1),perm(swap(:,1))')) ...
                - Aerror(sub2ind(size(Aerror),swap(:,2),perm(swap(:,2))')));
            
            if m >= 0 break; end
            
            temp = perm(swap(ind,1));
            perm(swap(ind,1)) = perm(swap(ind,2));
            perm(swap(ind,2)) = temp;
        end
        
        %perm
        %sum(Aerror(sub2ind(size(Aerror),1:L,perm)))
        
        if sum(Aerror(sub2ind(size(Aerror),1:L,perm))) < CE
            CE = sum(Aerror(sub2ind(size(Aerror),1:L,perm)));
        end
    end
    
    CE = CE / N;
    
end

end