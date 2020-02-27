function [atom_bag_list, N_added, DG0, Uo, Sigo, Vo] = bag_1d_svd(grad, x, N_bag, S0, tau, norm_type, w, on_boundary, US, SigS, VS)
%BAG_1d_SVD bagging with vectors using subspace bagging
%
%Input:
%   grad:       Current gradient at point x
%   x:          Position x, optimal with respect to its atoms from last BASH
%   N_bag:      Number of desired atoms in bag
%   S0:         Active set of atoms found in previous BASH
%   tau:        Radius of Atomic Ball
%   norm_type:  specific norm under discussion
%   w:          Additional weights for the norm if necessary
%   on_boundary: whether there is saturation in any coordinate
%   US SigS VS': Thin SVD of [S]
%Output:
%   atom_bag_list:  Atoms in the bag to be added to active set
%   DG0:            Duality gap of the best atom
%   N_added:        Number of atoms that successfully entered the bag
%   Uo*Sigo*Vo':    Thin SVD of [S0 S_bag], augmented atomic set
atom_bag_list = [];
DG_min = 1e-5;
DG_angle_min = 1e-9;


%number of atoms previously in bag (with independent directions)
N_prev = size(S0, 2) - on_boundary;

r = size(S0, 2);

%[US, SigS, VS] = svd(S, 'econ');
%W = tau*S - x;
%recenter by -x

%w is used for multiple things.
%w is the weights
%w_dir is the directin tau*s - x
%W is the matrix of w_dir's, which is not explicitly formed
if isempty(S0)
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
        last_atom = S0(:, r);
        s_aug = ones(1, r);
        s_aug(r) = 0;
        [US, SigS, VS] = atom_svd_bash(S0, US, SigS, VS, s_aug, 0);
        b = ones(r-1, 1);
    else
        b = ones(r, 1);
    end
    a = -x;
    [U, Sig, V] = rank_one_svd_update(US, tau*SigS, VS, a, b, 0);        
end

G0 = -grad;
G = G0;
S_total = S0;

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
    %LMO: find the best atom
    %a = sign(G);    %best atom
    a = LMO_1d(G, norm_type, w);
    w_dir = tau*a - x;      %actual direction of travel
    DG = real(G0'*w_dir);   %was using G as a reference. That's a serious bug.
    DG_angle = DG / (norm(G0) * norm(w_dir));
    
    if N_added == 0
        DG0 = DG;
    end
    
    %Update the thin SVD describing the atoms
    %Repeat bagging process until convergence occurs (empty or full bag)
    if (DG > DG_min) && (DG_angle > DG_angle_min)
        if isempty(U)
            Sigp = norm(w_dir);
            Up = w_dir / Sigp;
            Vp = 1;
        else
            %add a new column
            if mod(N_added, 10) == 9
                [Up, Sigp, Vp] = add_svd_column(U, Sig, V, w_dir, 1);
            else
                [Up, Sigp, Vp] = add_svd_column(U, Sig, V, w_dir, 0);
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
if ~isempty(S0)
    [Uo, Sigo, Vo] = rank_one_svd_update(U, Sig, V, x, b, 0);
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
    %add back the last atom
    [Uo, Sigo, Vo] = add_svd_column(Uo, Sigo, Vo, last_atom, 0, N_prev);    
end

end

