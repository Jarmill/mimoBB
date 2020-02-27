function [Up, Sigp, Vp] = atom_svd_bash(S, U, Sig, V, s_aug, force_orth)
%atom_svd_bash updates the the thin svd of the set of atoms U*S*V'
%modifications are in the array s, filled with signs. Used after the bash
%phase, when an optimal loading over atoms are attained.
%also modifies associated properties S, AS, and K_int
%
%Input:
%   S: current set of atoms
%   U*S*V': Thin svd of active set of atoms + bagged atoms
%   s_aug:  List of signs that reflect changes in atoms. augmented by bag
%               s = 1   Atom is kept, no change
%               s = 0   Atom is dropped, delete column
%               s = -1  Atom flips sign (if A is centrally symmetric)
%   force_orth: force orthogonality in SVD update
%dimensions of S (influence)
[d, r] = size(S);

%figure out what atoms get flipped and dropped

ind_flip = find(s_aug == -1);
num_flip = length(ind_flip);

ind_drop = find(s_aug == 0);
num_drop = length(ind_drop);
%need to also figure out updates of S, AS, and K_int

if num_flip > 0
    %form the block matrices for flipping
    A_flip = -2*S(:, ind_flip);
    i_flip = ind_flip;    
    j_flip = 1:num_flip;
    v_flip = ones(size(i_flip));
    B_flip = sparse(i_flip, j_flip, v_flip, r, num_flip);

    [U, Sig, V] = svd_update(U, Sig, V, A_flip, B_flip, force_orth);
end

%now the block matrices for deleting
if num_drop > 0
    A_drop = -S(:, ind_drop);
    i_drop = ind_drop;
    j_drop = 1:num_drop;
    v_drop = ones(size(i_drop));
    B_drop = sparse(i_drop, j_drop, v_drop, r, num_drop);

    [U, Sig, V] = svd_update(U, Sig, V, A_drop, B_drop, force_orth);
    %very delicate indexing over here.
    Up   = U(:, 1:(end-num_drop));
    Sigp = Sig(1:(end-num_drop), 1:(end-num_drop));
    Vp   = V(:, 1:(end-num_drop));
    Vp(ind_drop, :) = [];
else
    %no dropping necessary
    Up = U;
    Sigp = Sig;
    Vp = V;
end



