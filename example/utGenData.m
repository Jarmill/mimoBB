function [z,zn,noise,ny,nu,nx] = utGenData(sys,SNR,InputBWFraction,Ns)
% Experimental data realization.
% Generate data realization for a given system (sys), using an input that
% is WGN of power limited to InputBWFraction times the Nyquist frequency.
% Output noise of variance determined by SNR is added to the model outputs.

% poor excitation if InputBWFraction < 0.2;

[ny,nu] = size(sys);
nx = order(sys);
s = 10.^(SNR/10);
u = idfilt(randn(Ns,nu),[0 pi*InputBWFraction]);
y = lsim(sys,u); % note: zero IC
n = randn(Ns,ny);
noise = zeros(size(y));
for ky = 1:ny
   noise(:,ky) = sqrt(var(y(:,ky))/s)*n(:,ky);
end
if ny==2
   noise = noise*[1 0.2; 0.2 0.5];
end

yn = y + noise;
z = iddata(y,u,1);
zn = iddata(yn,u,1);
