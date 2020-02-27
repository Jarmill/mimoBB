function [h,p,scales, groups, f] = createAtoms(NumSamples,opt,G)
%% Create Rational Atoms
% Impulse response matrix, where each column corresponds to N-length
% impulse response of:
% +/-1/(z-p), if p is real
% +/-(1/(z-p) + 1/(z-p')), +/-(-j/(z-p) + j/(z-p')), is p is complex
% Impulse responses are scaled to unity nuclear norm.

% r1, r2: radius limits 
% MaxAngleFraction: (0,1), +ve portion of unit cirle to include centered
% around positive real axis. Example: 0.25 means angles between +/-(pi/2).

% Let p0 be vector of randomly generated atoms ; length(p0) = NumAtoms
% length(p) is 2*numel(real_p0) + 4*numel(complex_p0). This is to enforce
% symmetry.
%
%
% Output:
%   h:      impulse response array (time)
%   p:      set of poles
%   scales: scales on poles to ensure unit hankel norm
%   groups: groups of associated (complex conjugate) poles
%   f:      frequency response array (freq)


groups = [];
if opt.Type=="TC"
   h = [];
   for ct = 1:numel(opt.Alpha)
      h = [h,kernelSS(opt.Alpha(ct),NumSamples)];
   end
   scales = [];
   p = 1:size(h,2);
   return
elseif any(opt.Type==["spline","random"])
   h = createSSAtoms(NumSamples,opt);
   scales = [];
   p = 1:size(h,2);   
   return;
end

C = opt.IncludeConstant;
p = uniform_over_ring_sector(opt.r1,opt.r2,opt.NumAtoms,opt.phi1,opt.phi2,C);
if C
   pc = p(end);
   p = p(1:end-1);
end
p = p(:).'; % row vector
Ir = imag(p)==0;
nr0 = sum(Ir); nc0 = numel(p)-nr0;
%[scales,L] = getScales(p,NumSamples); % nr+2*nc elements
%p2 = [p(Ir),-p(Ir),p(~Ir),-p(~Ir)];
% duplicate scales for negative atoms.
[scales,L] = getScales(p,NumSamples); % nr0+2*nc0 elements
p_real = [p(Ir),-p(Ir)];
p_comp_0 = [p(~Ir),-p(~Ir)];
p_comp = [p(~Ir),conj(p(~Ir)),-p(~Ir),conj(-p(~Ir))];
p2 = [p_real, p_comp_0];
p = [p_real, p_comp];

% duplicate scales for negative atoms.
scales = [scales(1:nr0),scales(1:nr0),scales(nr0+1:end),scales(nr0+1:end)];

%Time response
k = numel(p2);
h = [zeros(1,k);ones(1,k);cumprod(p2(ones(1,NumSamples-2),:))];
Ic = (2*nr0+1):numel(p2);
hc1 = 2*real(h(:,Ic));
hc1p = hc1(:,1:nc0);
hc1n = hc1(:,nc0+1:end);
hc2 = -2*imag(h(:,Ic));
hc2p = hc2(:,1:nc0);
hc2n = hc2(:,nc0+1:end);
h = [h(:,1:2*nr0),hc1p,hc2p,hc1n,hc2n];
h = h.*scales;

%Frequency response
w = G.Frequency;
NumSamplesF = numel(w);
z = ltipack.utGetComplexFrequencies(w,1);
h1 = 1./(z-p(Ir));
h1a = 1./(z+p(Ir));
Dr = z.^2-2*real(p(~Ir)).*z+abs(p(~Ir)).^2;
Dra = z.^2+2*real(p(~Ir)).*z+abs(p(~Ir)).^2;
h2 = 2*imag(p(~Ir))./Dr;
h2a = 2*imag(p(~Ir))./Dra;
h3 = 2*(z-real(p(~Ir)))./Dr;
h3a = 2*(z+real(p(~Ir)))./Dra;
f = [h1,h1a,h3,h2,h3a,h2a];
f = f.*scales;

groups = [1:nr0, nr0+(1:nr0), 2*nr0 + kron([1,1], 1:nc0), 2*nr0 + nc0 + kron([1,1], 1:nc0)];

if C
   p = [p, pc];
   h = [h, ones(size(h,1),1)/L];
   groups = [groups, 2*(nr0+nc0)+1];
   f = [f, getScales(1,NumSamplesF)*1./(z-1)];
end

return

%alternatively use frequency coordinates available from sample

%I have no idea if this is correct, really hope it is.
%Should we be using fftshift?

w_base = exp(-2*1j*pi/NumSamples);

%I still have no idea
%w = linspace(0, NumSamples-1, NumSamples)'/NumSamples;
w = w_base.^((0:(NumSamples-1))');
ew = exp(-2*pi*1.0j*w);

f_real_denom = 1 - ew*p_real;
f_comp_denom = 1 - 2*ew*real(exp(-p_comp_0)) + ew*real(exp(-2*p_comp_0));
f_sin_num    = ew*imag(exp(-p_comp_0));
f_cos_num    = 1 - ew*real(exp(-p_comp_0));

f_real = 1./f_real_denom;
f_sin = f_sin_num./f_comp_denom;
f_cos = f_cos_num./f_comp_denom;

f = [f_real, f_cos(:, 1:nc0), f_sin(:, 1:nc0), f_cos(:,  nc0+1:end), f_sin(:, nc0+1:end)];
f = f.*scales;
%which atoms are associated together?
%complex poles will have their complex conjugate present
groups = [1:nr0, nr0+(1:nr0), 2*nr0 + kron([1,1], 1:nc0), 2*nr0 + nc0 + kron([1,1], 1:nc0)];

if C
   p = [p, pc];
   h = [h, ones(size(h,1),1)/L];
   groups = [groups, 2*(nr0+nc0)+1];
   f = [f, [L; zeros(size(f,1)-1,1)]];
end

% delays
hh = [diag(ones(1,10),-1);zeros(NumSamples-11,11)];
hh = hh(:,1:end-1)./(1:10);
%h = [h,hh];
%p = [p, zeros(1,10)];


