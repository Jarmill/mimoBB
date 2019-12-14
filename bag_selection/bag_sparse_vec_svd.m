function [atom_bag_list, N_added] = bag_sparse_vec_svd(grad, x, on_boundary, N_bag, CS, w)
%BAG_SPARSE_VEC_SVD performs subspace bagging over the set of sparse
%vectors. This includes the weighted simplex/L1-ball. Does not actually
%perform an svd, instead uses the structure of the simplex/L1-ball and easy
%formulas for vector rejection on these faces.
%
%Input:
%   grad:       gradient at current x
%   x:          Current vector x = S*c, regressor of interest
%   on_boundary:If x has full atomic norm
%   N_bag:      number of elements that would like to be bagged
%   norm_type:  which norm or gauge to use as the atomic set
%   w:          Set of weights on reweighted heuristic, diagonal scaling
%
%Output:
%   atom_bag_list:  Set of atoms in the bag, to be added in BASH
%   N_added:        Number of atoms in the bag

%start by accounting for weighting
if nargin < 4 || isequal(w, 1)
    weighted = 0;
    grad_weighted = grad;
else
    weighted = 1;
    grad_weighted = grad ./ w;
end

%set up iterative subspace bagging and list of signs s. s is also the
%normal vector for the L1-ball or main face of the simplex (not sidewall)
G0 = -grad_weighted;
G = G0;
s = sign(x);
d = size(x, 1);

if on_boundary
    %start by forming the sparse vector
    i = [];
    j = [];
    v = [];
    N_added = 0;
    
    %main loop, keep adding dimensions until the bag is full
    while N_added < N_bag
        %find best index and respective value
        if CS
            [~, i_curr] = max(abs(G));
        else
            [~, i_curr] = max(G);
        end
        
        v_curr = sign(G(i_curr));
        s(i_curr) = v_curr;
        if weighted
            v_curr = v_curr / w(i_curr);
        end
        
        %if a valid choice was made
        if s(i_curr) == 0            
            %append to sparse matrix describing atoms
            N_added = N_added + 1;
            i = [i i_curr];
            j = [j N_added];
            v = [v v_curr];
            
            %update the 'negative gradient' based on the vector rejection
            %amount outside the subspace
            %this may be exactly equivalent to cluster bagging. I think it
            %is exactly equivalent to cluster bagging, as the L1 ball is
            %orthogonal and nice and pleasant compared to other atomic
            %domains. The orthogonality of the atomic set actually fuels 
            G(s ~= 0) = 0;
        else
            break
        end
            
    end
else
    %on interior of region, things are easy
    if CS
        %L1-ball
        [~, valid_dim] = maxk(abs(G), N_bag);
    else
        [~, valid_dim] = maxk(G, N_bag);
    end
    
    %only dimensions which are currently zero are eligible to be added
    added_dim = valid_dim(x(valid_dim) == 0);
    N_added = length(added_dim);
    
    %pull together the sparse vector
    i = added_dim;
    j = 1:N_added;
    v = sign(G(i));
    if weighted
        v = v ./ w(added_dim);
    end
end

atom_bag_list = sparse(i, j, v, d, N_added);
end