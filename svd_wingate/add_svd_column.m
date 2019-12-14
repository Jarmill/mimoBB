function [Ua, Siga, Va] = add_svd_column(U, S, V, c, force_orth, ind)
%add_svd_column add a new column to the thin svd of A = U S V'
%
%Input:
%   U S V':     SVD of matrix
%   c:          inserted column
%   force_orth: force orthogonality of the thin svd decomposition
%   ind:        What column to add the index. By default is the end.
%
%Output
%   Ua Siga Va': SVD of low rank update with column c added at index ind
d = size(V, 2);


b = zeros(d+1, 1);
z = zeros(1, d);

if nargin < 6
    %add last column
    ind = d+1;
    b(ind) = 1;
    Vz = [V; z];
else
    %add to middle
    if ind == d+1
        b(ind) = 1;
        Vz = [V; z];
    else
        b(ind+1) = 1;
        Vz = [V(1:ind, :); z; V((ind+1):end, :)]; 
    end
end    

%b(ind) = 1;


a = c;
[Ua, Siga, Va] = rank_one_svd_update(U, S, Vz, a, b, force_orth);