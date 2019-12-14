function A = broken_arrowhead(s, m, p)
%BROKEN_ARROWHEAD produces a sparse broken arrowhead matrix.
%
%Input:
%   s:  main diagonal of arrowhead (singular values)
%   m:  right row on the end
%   p:  entry at bottom right corner
%
%Output:
%   A:  broken arrowhead output matrix

%rank (size of output matrix)
r = length(s) + 1;

i_arrow = [1:r 1:(r-1)];
j_arrow = [1:r r*ones(1, r-1)];
v_arrow = [s; p; m];

A = sparse(i_arrow, j_arrow, v_arrow, r, r);
%A = spdiag([s p]', 0, r, r);
