function [atom_bag_list, N_added] = bag_l_inf(grad, x, N_bag, S, tau, on_boundary)
%BAG_L_INF bagging over the L_infinity ball (max ball/cube)
%dual norm is the l1 norm
%
%Input:
%   grad:   Current gradient at point x
%   x:      Position x, optimal with respect to its atoms from last BASH
%   N_bag:  Number of desired atoms in bag
%   S:      Active set of atoms found in previous BASH
%   tau:    Radius of L-Infinity Ball
%   on_boundary: whether there is saturation in any coordinate
%Output:
%   atom_bag_list:  Atoms in the bag to be added to active set
%   N_added:        Number of atoms that successfully entered the bag
atom_bag_list = [];
DG_min = 1e-4;

%number of atoms previously in bag (with independent directions)
N_prev = size(S, 2) - on_boundary;

if isempty(S)
    W = [];
else
    if on_boundary
        %the output of bash on the boundary is normal to all atoms, and
        %there is a rank drop due to the boundary constraint. Get rid of
        %the last column to compensate
        W = tau*S(:, 1:N_prev) - x;
    else
        W = tau*S - x;
    end
end

[Q, R] = qr(W);

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



%there's still lockout when i try and add more than one dimension. Fix this
%later.

while ~done
    a = sign(G);    %best atom
    w = tau*a - x;      %actual direction of travel
    DG = G0'*w;
    

    
    %really need a more efficient method than ismember
    %if G0'*(tau*a-x) > 1e-6 && isempty(S_total) || ~ismember(a', S_total', 'rows')
    if DG > DG_min
        if N_bag > 1
            %sequence of directions travelled
            %W = [W w];
            %should probably figure out SVD update for adding columns
            %[U, ~, ~] = svd(W, 'econ');
            %U(abs(U) < 1e-9) = 0;
            %G_proj = U*(U'*G0);
            %G = G - G_proj;
            %[Q, R] = qr(W, 0);
            %[Q, R] = qr(W);            
            [Q, R] = qrinsert(Q, R, N_prev + N_added + 1, w);
            Q_active = Q(:, 1:(N_prev + N_added + 1));
            %G_proj = Q*(Q'*G0);
            G_proj = Q_active*(Q_active'*G0);
            G = G0 - G_proj;
        end

        
        %there could be repeats. How to deal with repeats?

        atom_bag_list = [atom_bag_list a];
        %S_total = [S_total a];
        N_added = N_added+1;
    else
       done = 1;
    end
    
    if N_added >= N_bag
        done = 1;
    end
end



end

