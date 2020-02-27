function [atom_bag_list, N_added] = bag_complex_1d(grad, x, N_bag, tau, w)
%BAG_COMPLEX_1D complex numbers on L1 constraint
%Is a sum of L2 norms
%x is optimal over its nonzero dimensions from a previous iteration. Thus
%the use of theta (the lagrange multiplier).
%
%Input:
%   grad:   Gradient at current x
%   x:      Current point x
%   N_bag:  Number of new atoms to add if possible
%   tau:    Atomic norm constraint
%   w:      Weights on atoms from reweighted heuristic (optional)
%
%Output:
%   atom_bag_list:  Atoms in bag to be added to system, and tried in BASH
%   N_added:        Number of atoms successfully added to bag
if nargin < 5 || isequal(w, 1)
    weighted = 0;
    grad_weighted = grad;
else
    weighted = 1;
    grad_weighted = grad ./ w;
end

%find best atoms
absgrad = abs(grad_weighted);
    
[~, added_dim] = maxk(absgrad, N_bag);

%form the atoms into a matrix
i = added_dim;
N_added = length(added_dim);
j = 1:N_added;

if weighted
    w_curr = w(added_dim);
else
    w_curr = 1;
end
v = -grad(i) ./ abs(grad(i)) ./ w_curr;

%atoms that are added to the system
atom_bag_list_screen = sparse(i, j, v, size(x, 1), N_added);
dir_travel = tau*atom_bag_list_screen - x;
DG = real(dir_travel'*(-grad));
DG_best = max(DG);
%some real ugly code down here.
survive_index = find(DG > 1e-5);
N_added = length(survive_index);
atom_bag_list = atom_bag_list_screen(:, survive_index);
%DG = dir_travel'*grad;
%check duality gap

% else
%     [~, added_dim] = maxk(abs(grad_weighted(missing_coord_candidates)), N_bag);
%     
%     N_added = length(added_dim);
%     i = missing_coord_candidates(added_dim);
%     j = 1:N_added;
%     %v = ones(size(j));
%     %scale atoms in case of reweighted heuristic
%     if weighted
%         w_curr = w(added_dim);
%     else
%         w_curr = 1;
%     end
%     
%     if isreal(grad_weighted)
%         v = sign(-grad_weighted(i)) ./ w_curr;
%     else
%         v = -grad_weighted(i) ./ abs(grad_weighted(i)) ./ w_curr;
%     end
%     %atoms that are added to the system
%     atom_bag_list = sparse(i, j, v, size(x, 1), N_added);
%     
% end

end

