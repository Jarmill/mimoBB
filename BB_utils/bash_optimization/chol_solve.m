function [y] = chol_solve(L, r)
%CHOL_SOLVE Solves the system K y = r where K = L'*L (cholesky
%factorization). Is a very common task
%
%Input:
%   L:  Cholesky decomposition of the kernel L = chol(K)
%   r:  Right hand side of the expression
%
%Output:
%   y:  answer of the SPD system

Ly = L' \ r;
y  = L \ Ly;

end

