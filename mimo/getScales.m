function [scales,L] = getScales(p,N)
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

Ir = imag(p)==0;
scales1 = (1-p2(Ir))./(1-p2N(Ir).*p2(Ir));
scales1(I(Ir)) = 1/L;

a1_sq = real(a1.^2);
c1_sq = c1.^2;
s2 = 2*sqrt(2*Gam.*(abs(a1).^2-c1.^2));
scales2 = 1./sqrt(2*(a1_sq + c1_sq) + s2);
scales3 = 1./sqrt(2*(-a1_sq + c1_sq) + s2);
scales = [scales1,scales2(~Ir),scales3(~Ir)];

end
