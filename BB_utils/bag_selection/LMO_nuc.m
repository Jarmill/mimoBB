function [u, v] = LMO_nuc(G)
%LMO_NUC Linear Minimization Oracle for the Nuclear Norm
%
%Input:
%   G:  negative gradient or matrix reference

[u, s, v] = svds(G, 1);

end

