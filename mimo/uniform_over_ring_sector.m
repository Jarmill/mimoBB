% Generate n uniformly randomly distributed complex numbers inside the ring
% defined by rho1 <= |p| <= rho2 and sector of MinAngle to MaxAngle
% This code does not force set symmetry.


function p = uniform_over_ring_sector(rho1,rho2,n,MinAngle,MaxAngle,IncludeConstant)

Offset = nargin>5 && IncludeConstant;
n0 = n;
n = 10*n0;
%Range = (rand(1,n)-0.5);
Range = 0.5*rand(1,n); %conjugates will be sampled anyways
theta = MinAngle + (MaxAngle-MinAngle)*Range; 


%number of real poles
nr = floor(n/6);

%some poles should be purely real
if MaxAngle == 2*pi
    %if negative poles are allowed by sector bound,
    %then some real poles should be negative
    signflip = randi(2, 1, nr)-1;
    signflip(signflip == 1) = pi;
    theta(1:nr)  = signflip;
else
    theta(1:nr) = 0;
end

theta = theta(randperm(n)); % scramble 0 angle entries
r = sqrt(rand(1,n))*(rho2-rho1);
p = (r+rho1).*cos(theta)+1j*(r+rho1).*sin(theta);
I = abs(angle(p))<=MaxAngle & abs(angle(p))>=MinAngle & abs(p)<=rho2 & abs(p)>=rho1;
%dd = -cos(angle(log(p))); I = I & dd<0.2;
p = p(I);
p = p(1:min(n0,end));
if Offset
   p(end) = 1;
end

%{
figure
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
