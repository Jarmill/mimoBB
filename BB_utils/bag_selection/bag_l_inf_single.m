function [atom_bag_list, N_added] = bag_l_inf_single(grad, x, N_bag, S, tau)
%BAG_L_INF bagging over the L_infinity ball (max ball/cube)
%dual norm is the l1 norm
%
%Input:
%   grad:   Current gradient at point x
%   x:      Position x, optimal with respect to its atoms from last BASH
%   N_bag:  Number of desired atoms in bag
%   S:      Active set of atoms found in previous BASH
%   tau:    Radius of L-Infinity Ball
%Output:
%   atom_bag_list:  Atoms in the bag to be added to active set
%   N_added:        Number of atoms that successfully entered the bag

atom_bag_list = [];

G0 = -grad;
G = G0;
%S_total = S;

N_added = 0;
done = 0;


%there's still lockout when i try and add more than one dimension. Fix this
%later.
a = sign(G);
DG = -grad'*(tau*a-x);
if DG > 1e-4
    atom_bag_list = [atom_bag_list a];
    %S_total = [S_total a];
    N_added = 1;
end


%there could be repeats. How to deal with repeats?

end

