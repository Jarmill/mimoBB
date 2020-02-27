% Script to find the scales for atoms when atoms are expressed solely by
% their impulse response. For each atom, the scale normalizes the
% nuclear norm of Hp, where Hp is the largest square Hankel matrix one can
% construct using given response.

function scale_r = getSSScales(h)
% h: impulse response matrix

[N_hor,NumAtoms1,NumAtoms2] = size(h);
if mod(N_hor,2) == 0
   N_odd =  2*ceil(N_hor/2)-1;   % Nearest odd number <= N, so that Hankel is square
   L = (N_odd+1)/2;              % Dimension of the square hankel matrix
else
   N_odd =  2*ceil(N_hor/2)-1;   % Nearest odd number <= N, so that Hankel is square
   L = (N_odd+1)/2-1;            % Dimension of the square hankel matrix
end

scale_r = zeros(NumAtoms1,NumAtoms2);
for k = 1:NumAtoms1*NumAtoms2
   H = hankel(h(2:L,k),h(L:end,k));
   scale_r(k) = 1./norm_nuc(H);
end

