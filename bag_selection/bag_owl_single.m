function [atom_bag_list, N_added] = bag_owl(grad, x, w, N_bag, S, tau)
%BAG_OWL bagging over OWL norm (ordered weighted L1 norm)
%dual norm is the I have no idea.
%Input:
%   grad:   Current gradient at point x
%   x:      Position x, optimal with respect to its atoms from last BASH
%   w:      Weights on dimensions for owl norm (weakly decreasing sequence)
%   N_bag:  Number of desired atoms in bag
%   S:      Active set of atoms found in previous BASH
%   tau:    Radius of OWL-ball
%
%Output:
%   atom_bag_list:  Atoms in the bag to be added to active set
%   N_added:        Number of atoms that successfully entered the bag
atom_bag_list = [];
% if isempty(S)
%     W = [];
% else
%     W = tau*S - x;
% end

G0 = -grad;
G = G0;

%iterative gram-schmidt type of scheme. Orthogonal Bagging
%will add Cluster bagging later
N_added = 0;

ws = cumsum(w);         %weights are a weakly decreasing sequence
wr = 1./ws;             %reciprocol of the cumulative sum, is the weights 
                        %on each atom for a set of coordinates chosen


%find the best atom in the set
    
%find the indices with the highest gradient
%[~, i_max] = maxk(abs(G), N_weights);
[G_abs_sort, i_max] = sort(abs(G), 'descend');

%G_abs_sort = abs(G(i_max));
%get the cumulative sums of the greatest absolute values
G_abs_sum = cumsum(G_abs_sort);

%compare the weighted sum of using k components. The more components
%used, the higher the cumulative sum, but the lower the associated
%atomic weight (sum of reciprocals of weights w to get to that point)
%
%Use the linear minimization oracle to decide (subgradient of dual OWL
%norm)
LMO = wr .* G_abs_sum;
[~, N_components] = max(LMO); %how many components to add
ind_added = i_max(1:N_components); %the indices corresponding to added components


%the final atom is the signs of G's components weighted by the atomic
%weight of adding that many components. Kind of confusing, but easier
%than the OWL paper's hadamard product and permutation and orderings.
%a = sparse(length(x), 1);
a = zeros(size(x));
a(ind_added) = sign(G(ind_added)) * wr(N_components);

                        
w_dir = tau*a - x;                        
DG = G0'*w_dir;
if DG > 1e-4
    atom_bag_list = [atom_bag_list a];
    %S_total = [S_total a];
%else
%    done = 1;
end

end

