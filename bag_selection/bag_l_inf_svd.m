function [atom_bag_list, N_added, Uo, Sigo, Vo] = bag_l_inf_svd(grad, x, N_bag, S, tau, on_boundary, US, SigS, VS)
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
%   Uo*Sigo*Vo':    Thin SVD of [S0 S_bag], augmented atomic set
atom_bag_list = [];
DG_min = 1e-4;

%number of atoms previously in bag (with independent directions)
N_prev = size(S, 2) - on_boundary;

r = size(S, 2);

%[US, SigS, VS] = svd(S, 'econ');
%W = tau*S - x;
%recenter by -x
if isempty(S)
    U = [];
    Sig = [];
    V = [];
else
    %when x is on the boundary, then it is optimal with respect to its
    %atoms, and the gradient points outwards in line with the face normal.
    %Linear independence among directions of motion (travel towards
    %tau-scaled atoms) is lost. This messes with projections and subspaces.
    %If x is on the boundary, remove the last atom from S, and add it in at
    %the end.
    if on_boundary
        last_atom = S(:, r);
        s_aug = ones(1, r);
        s_aug(r) = 0;
        [US, SigS, VS] = atom_svd_bash(S, US, SigS, VS, s_aug, 0);
        b = ones(r-1, 1);
    else
        b = ones(r, 1);
    end
    a = -x;
    [U, Sig, V] = rank_one_svd_update(US, tau*SigS, VS, a, b, 0);        
end

%[Q, R] = qr(W);
%[U, Sig, V] = svd(W, 'econ');

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
    DG = G'*w;
    

    
    %really need a more efficient method than ismember
    %if G0'*(tau*a-x) > 1e-6 && isempty(S_total) || ~ismember(a', S_total', 'rows')
    if DG > DG_min
        if isempty(U)
            Sigp = norm(w);
            Up = w / Sigp;
            Vp = 1;
        else
            %add a new column
            if mod(N_added, 10) == 9
                [Up, Sigp, Vp] = add_svd_column(U, Sig, V, w, 1);
            else
                [Up, Sigp, Vp] = add_svd_column(U, Sig, V, w, 0);
            end
        end
        G_proj = Up*(Up'*G0);
        G = G0 - G_proj;

        U = Up;
        Sig = Sigp;
        V = Vp;    
        
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


b = ones(size(V, 1), 1);
if ~isempty(S)
    [Uo, Sigo, Vo] = rank_one_svd_update(U, Sig, V, x, b, 1);
else
    Uo = U;
    Sigo = Sig;
    Vo = V;
    
    %re-add last atom to system
end
Sigo = Sigo/tau;

%last column loses linear independence on the boundary, screwing over
%projections. Restore it and swap it to its rightful position. Ordering of
%atoms in the bag doesn't matter to bash, as it kills atoms by monitoring
%lockout. 
if on_boundary
    [Uo, Sigo, Vo] = add_svd_column(Uo, Sigo, Vo, last_atom, 0);
    %swap_columns(Vo, r, r+N_added); %actually need to swap rows, not cols.
    Vo([r, r+N_added], :) = Vo([r+N_added, r], :);
    atom_bag_list = circshift(atom_bag_list, -1, 2);
end

end

