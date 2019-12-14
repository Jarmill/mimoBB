% Script to find the scales for atoms. 

% Author : Burak Yilmaz 
% Last Update :  10/27/2014

% Important : N is the horizon length, not the dimenstion of the Hankel 

% The code is little different than what is derived in the Paper Appendix.
% But the result is the same. For each atom, the scale normalizes the
% nuclear norm of Hp, where Hp is the largest square Hankel matrix one can
% construct using; the N_hor length impulse response of the atom after
% discarding the initial zero. 

function [scale_r,scale_c] = get_scales(p,N_hor)

    if mod(N_hor,2) == 0
    N_odd =  2*ceil(N_hor/2)-1;   % Nearest odd number <= N, so that Hankel is square
    L = (N_odd+1)/2;              % Dimension of the square hankel matrix
    else
    N_odd =  2*ceil(N_hor/2)-1;   % Nearest odd number <= N, so that Hankel is square
    L = (N_odd+1)/2-1;              % Dimension of the square hankel matrix
    end
    p2 = p.^2;
    p2N = p2.^L;
  
    pabs2 = abs(p).^2;
    pabs2N = pabs2.^L;
    a1 = (1-p2N)./(1-p2);
    c1 = (1-pabs2N)./(1-pabs2);
    M = (real(a1)-c1-real(p2.*conj(p2N).*a1)+pabs2.*pabs2N.*c1)./(1-pabs2.^2);

    prod_eig = 2*((abs(a1).^2-c1.^2).*M);

    a1_sq = real(a1.^2);
    c1_sq = c1.^2;

    sum_eig_r = 2*(a1_sq+c1_sq);
    sum_eig_c = 2*(c1_sq-a1_sq);
    norm_nuc_r = sqrt(sum_eig_r + 2*sqrt(prod_eig));
    norm_nuc_c = sqrt(sum_eig_c + 2*sqrt(prod_eig));
    scale_r = 1./norm_nuc_r;
    scale_c = 1./norm_nuc_c;

end