% Generate n uniformly randomly distributed complex numbers inside the ring
% defined by rho1 <= |p| <= rho2 and sector of +/- MaxAngle*2*pi


function p = uniform_over_ring_sector(rho1,rho2,n,MaxAngle)


u = (rho2 - rho1) * rand(n, 1) + rho1;
v = (rho1 + rho2) * rand(n, 1);

r = zeros(1, n);
r(v < u)  = u(v < u);
r(v >= u) = rho1 + rho2 - u(v >= u);

theta = 2*MaxAngle*(rand(1, n) - 0.5);

p = r .* (cos(theta) + 1j*sin(theta));

%Old code:
% n0 = n;
% n = 10*n0;
% Ang = 2*pi*MaxAngle;
% theta = randn(1,n)*(Ang);
% r = sqrt(randn(1,n))*(rho2-rho1);
% p = (r+rho1).*cos(theta)+1j*(r+rho1).*sin(theta);
% I = abs(angle(p))<=Ang & abs(p)<=rho2 & abs(p)>=rho1;
% p = p(I);
% p = p(1:min(n0,end));


%{
%figure
r = exp(1i*2*pi*(0:100)/100);
plot(real(p),imag(p),'b*',...
   real(r*rho1),imag(r*rho1),'k--',real(r*rho2),imag(r*rho2),'k--',...
   real(r),imag(r),'k')
%hold on
axis('square')
xlim([-1.1 1.1])
ylim([-1.1 1.1])
%}
end

