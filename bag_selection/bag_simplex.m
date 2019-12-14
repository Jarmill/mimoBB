function [atom_bag_list, N_added] = bag_simplex(grad, x, N_bag, w)
%BAG_SIMPLEX bagging over the simplex
%dual gauge is the minimum entry
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

%compensate for reweighted heuristic
if nargin < 4 || isequal(w, 1)
    weighted = 0;
    grad_weighted = grad;
else
    weighted = 1;
    grad_weighted = grad ./ w;
end

theta = -min(grad_weighted(x~=0));

if isempty(theta)
    missing_coord_candidates = find(grad_weighted < 0);
else
    missing_coord_candidates = find(grad_weighted < -theta);
end

if isempty(missing_coord_candidates)
    atom_bag_list = [];
    N_added = 0;
    DG = 0;
else
    [~, added_dim] = maxk(abs(grad_weighted(missing_coord_candidates)), N_bag);
    
    N_added = length(added_dim);
    i = missing_coord_candidates(added_dim);
    j = 1:N_added;
    %scale atoms in case of reweighted heuristic
    if weighted
        v = sign(-grad_weighted(i)) ./ w(added_dim);
    else
        v = sign(-grad_weighted(i));
    end

    %atoms that are added to the system
    atom_bag_list = sparse(i, j, v, size(x, 1), N_added);
end

end

