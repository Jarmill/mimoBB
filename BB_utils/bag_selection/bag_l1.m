function [atom_bag_list, N_added] = bag_l1(grad, x, N_bag, w)
%BAG_L1 bagging over the L1-ball
%dual norm is the l_inf (max)
%x is optimal over its nonzero dimensions from a previous iteration. Thus
%the use of theta (the lagrange multiplier).
%
%Input:
%   grad:   Gradient at current x
%   x:      Current point x
%   N_bag:  Number of new atoms to add if possible
%   w:      Weights on atoms from reweighted heuristic (optional)
%
%Output:
%   atom_bag_list:  Atoms in bag to be added to system, and tried in BASH
%   N_added:        Number of atoms successfully added to bag
if nargin < 4 || isequal(w, 1)
    weighted = 0;
    grad_weighted = grad;
else
    weighted = 1;
    grad_weighted = grad ./ w;
end


theta = max(abs(grad_weighted(x~=0)));
if isempty(theta)
    missing_coord_candidates = find(grad_weighted);
else
    missing_coord_candidates = find(abs(grad_weighted) > theta);
end

if isempty(missing_coord_candidates)
    atom_bag_list = [];
    N_added = 0;
else
    [~, added_dim] = maxk(abs(grad_weighted(missing_coord_candidates)), N_bag);
    
    N_added = length(added_dim);
    i = missing_coord_candidates(added_dim);
    j = 1:N_added;
    %v = ones(size(j));
    %scale atoms in case of reweighted heuristic
    if weighted
        w_curr = w(added_dim);
    else
        w_curr = 1;
    end
    
    if isreal(grad_weighted)
        v = sign(-grad_weighted(i)) ./ w_curr;
    else
        v = -grad_weighted(i) ./ abs(grad_weighted(i)) ./ w_curr;
    end
    %atoms that are added to the system
    atom_bag_list = sparse(i, j, v, size(x, 1), N_added);
    
end

end

