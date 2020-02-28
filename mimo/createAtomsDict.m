function [h,p, f, scales, groups] = createAtomsDict(p,w, C, Ns)
%CREATEATOMSDICT Summary of this function goes here
%   Detailed explanation goes here

%
%Input:
%   p:  Poles to generate responses
%   w:  Frequencies to sample in response
%   Ns: Number of samples (assume unit sample rate) 
%   C:  Include constant?
%
%Output:
%   h:  Impulse response
%   f:  Frequency response
%   scales: Hankel Norms of response (normalization factors)
%   groups: Which poles are clustered together?
if C
   pc = p(end);
   p = p(1:end-1);
end


p = p(:).'; % row vector
Ir = imag(p)==0;
nr0 = sum(Ir); nc0 = numel(p)-nr0;
%[scales,L] = getScales(p,Ns); % nr+2*nc elements
%p2 = [p(Ir),-p(Ir),p(~Ir),-p(~Ir)];
% duplicate scales for negative atoms.
[scales,L] = getScales(p,Ns); % nr0+2*nc0 elements
p_real = [p(Ir),-p(Ir)];

%p_comp_0 = [p(~Ir),-p(~Ir)];
%p_comp = [p(~Ir),conj(p(~Ir)),-p(~Ir),conj(-p(~Ir))];
p_comp_0 = [p(~Ir),-p(~Ir)];
p_comp = [ interleave2(p(~Ir),conj(p(~Ir))), interleave2(-p(~Ir),conj(-p(~Ir)))];

p2 = [p_real, p_comp_0];
p = [p_real, p_comp];

% duplicate scales for negative atoms.
scales = [scales(1:nr0),scales(1:nr0),scales(nr0+1:end),scales(nr0+1:end)];

%Time response
k = numel(p2);
h = [zeros(1,k);ones(1,k);cumprod(p2(ones(1,Ns-2),:))];
Ic = (2*nr0+1):numel(p2);
hc1 = 2*real(h(:,Ic));
hc1p = hc1(:,1:nc0);
hc1n = hc1(:,nc0+1:end);
hc2 = -2*imag(h(:,Ic));
hc2p = hc2(:,1:nc0);
hc2n = hc2(:,nc0+1:end);


h = [h(:,1:2*nr0),interleave2(hc1p,hc2p, 'col'),interleave2(hc1n,hc2n, 'col')];
%h = [h(:,1:2*nr0),hc1p,hc2p,hc1n,hc2n];
h = h.*scales;

%Frequency response
NsF = numel(w);
z = ltipack.utGetComplexFrequencies(w,1);
h1 = 1./(z-p(Ir));
h1a = 1./(z+p(Ir));
Dr = z.^2-2*real(p(~Ir)).*z+abs(p(~Ir)).^2;
Dra = z.^2+2*real(p(~Ir)).*z+abs(p(~Ir)).^2;
h2 = 2*imag(p(~Ir))./Dr;
h2a = 2*imag(p(~Ir))./Dra;
h3 = 2*(z-real(p(~Ir)))./Dr;
h3a = 2*(z+real(p(~Ir)))./Dra;
%f = [h1,h1a,h3,h2,h3a,h2a];
f = [h1, h1a, interleave2(h3,h2, 'col'), interleave2(h3a, h2a, 'col')];
f = f.*scales;

%groups = [1:nr0, nr0+(1:nr0), 2*nr0 + kron([1,1], 1:nc0), 2*nr0 + nc0 + kron([1,1], 1:nc0)];
groups = [1:nr0, nr0+(1:nr0), 2*nr0 + kron(1:nc0, [1,1]), 2*nr0 + nc0 + kron(1:nc0, [1,1])];

if C
   p = [p, pc];
   h = [h, ones(size(h,1),1)/L];
   groups = [groups, 2*(nr0+nc0)+1];
   f = [f, getScales(1,NsF)*1./(z-1)];
end

end

