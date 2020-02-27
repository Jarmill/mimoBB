function [atom_bag_list, N_added] = bag_l_p_single(grad, x, N_bag, S, tau, p)
%BAG_L_P_SINGLE bagging over the L_p ball
%dual norm is the L_q norm where 1/p + 1/q =1
%adding one atom at a time
%valid for p >= 1
%
%Input:
%   grad:   Current gradient at point x
%   x:      Position x, optimal with respect to its atoms from last BASH
%   N_bag:  Number of desired atoms in bag
%   S:      Active set of atoms found in previous BASH
%   tau:    Radius of L-Infinity Ball
%   p:      value of p for L_p norm
%Output:
%   atom_bag_list:  Atoms in the bag to be added to active set
%   N_added:        Number of atoms that successfully entered the bag
atom_bag_list = [];
DG_min = 1e-4;

if isempty(S)
    W = [];
else
    W = tau*S - x;
end

G0 = -grad;
G = G0;
S_total = S;

%iterative gram-schmidt type of scheme. Orthogonal Bagging
%will add Cluster bagging later
N_added = 0;
done = 0;

%add in the best atom not currently in the active set
%if an atom in the set is attempted to be added, terminate (or add a
%suboptimal atom? maybe that's also a bad idea).

%Figure out best atom by LMO
sg = sign(G);
ag = abs(G);

num_a = ag.^(1/(p-1));
den_a = norm(num_a, p);

a = sg .* num_a / den_a;
w = tau*a - x;
DG = G'*w;

if DG > 1e-4
    atom_bag_list = [atom_bag_list a];
    S_total = [S_total a];
    N_added = N_added+1;
end

%there's still lockout when i try and add more than one dimension. Fix this
%later.

% while ~done
%     a = sign(G);    %best atom
%     w = tau*a - x;      %actual direction of travel
%     DG = G0'*w;
%     
% 
%     
%     %really need a more efficient method than ismember
%     %if G0'*(tau*a-x) > 1e-6 && isempty(S_total) || ~ismember(a', S_total', 'rows')
%     if DG > DG_min
%         if N_bag > 1
%             %sequence of directions travelled
%             W = [W w];
%             %should probably figure out SVD update for adding columns
%             [U, ~, ~] = svd(W, 'econ');
%             U(abs(U) < 1e-9) = 0;
%             G_proj = U*(U'*G0);
%             %G = G - G_proj;
%             G = G0 - G_proj;
%         end
% 
%         
%         %there could be repeats. How to deal with repeats?
% 
%         atom_bag_list = [atom_bag_list a];
%         %S_total = [S_total a];
%         N_added = N_added+1;
%     else
%        done = 1;
%     end
%     
%     if N_added >= N_bag
%         done = 1;
%     end
% end



end

