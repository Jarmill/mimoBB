function [Ud, Sigd, Vd] = delete_svd_column(U, Sig, V, ind, c)
%delete_svd_column delete the ith column from the thin svd of A = U S V'
%
%Input:
%   U Sig V':   Thin SVD of matrix of atoms S
%   i:          Index of column to be dropped
%   c:          Entries column
%
%Output:
%
%   Ud Sigd Vd': Thin SVD of matrix of atoms S with ind's column removed
%really should find a better way to do this

r = size(V, 2);
%d = size(U, 1);


a = c;
b = zeros(r, 1);
b(ind) = -1;
[Ud, Sigd, Vd] = rank_one_svd_update(U, Sig, V, a, b, 0);

Ud = Ud(:, 1:(end-1));
Sigd = Sigd((1:end-1), 1:(end-1));
Vd = Vd(:, 1:(end-1));

Vd(ind, :) = [];