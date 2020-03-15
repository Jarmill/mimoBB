function [scales,L] = getScales2(p,N)
% Obtain scales such that nuclear norm of hankel(impulse(sys(p))) is 1.

if mod(N,2) == 0
   N_odd =  2*ceil(N/2)-1;   % Nearest odd number <= N, so that Hankel is square
   L = (N_odd+1)/2;              % Dimension of the square hankel matrix
else
   N_odd =  2*ceil(N/2)-1;   % Nearest odd number <= N, so that Hankel is square
   L = (N_odd+1)/2-1;              % Dimension of the square hankel matrix
end
p2 = p.^2;
p2N = p2.^L;

pabs2 = abs(p).^2;
pabs2N = pabs2.^L;
a1 = (1-p2N)./(1-p2);
c1 = (1-pabs2N)./(1-pabs2);
I = p2==1;
Gam = (real(a1)-c1-real(p2.*conj(p2N).*a1)+pabs2.*pabs2N.*c1)./(1-pabs2.^2);
Gam(I) = 0;

%Ir = imag(p)==0;
%scales1 = (1-p2(Ir))./(1-p2N(Ir).*p2(Ir));
scales = (1-pabs2)./(1-pabs2N.*pabs2);
scales(I) = 1/L;
