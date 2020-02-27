function G = grad_nuc(M, d, ind, delta)
%GRAD_NUC gradient of the nuclear norm. Can be made more efficient.
%   M:      Matrix M as input
%   d:      values  where data is known
%   ind:    Indices where data is known
%   delta:  L2 regularization paramter

proj_term = M(ind);
proj_term(ind) = proj_term(ind) - d;

G = proj_term + delta*M;


end