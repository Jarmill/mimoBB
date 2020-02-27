function [h,p,scales] = createAtomsFD(w,opt)
%% Create Rational Atoms
% Impulse response matrix, where each column corresponds to N-length
% impulse response of:
% +/-1/(z-p), if p is real
% +/-(1/(z-p) + 1/(z-p')), +/-(-j/(z-p) + j/(z-p')), is p is complex
% Scaling applied such that the corresponding impulse responses are unit
% nuclear norm in rational atom case.
%
% r1, r2: radius limits
% MaxAngleFraction: (0,1), +ve portion of unit cirle to include centered
% around positive real axis. Example: 0.25 means angles between +/-(pi/2).

% Let p0 be vector of randomly generated atoms ; length(p0) = NumAtoms
% length(p) is 2*numel(real_p0) + 4*numel(complex_p0). This is to enforce
% symmetry.

NumSamples = numel(w);
z = ltipack.utGetComplexFrequencies(w,1);
C = opt.IncludeConstant;

if any(opt.Type==["spline","random","TC"])
   opt.IncludeConstant = false;
   if opt.Type=="TC"
      h = [];
      for ct = 1:numel(opt.Alpha)
         K = kernelSS(opt.Alpha(ct),NumSamples);
         [K,~] = eig(K);
         h = [h, K];
      end
   else
      h = createSSAtoms(NumSamples,opt); h = h(:,1:end-10); % exclude delay elements
   end
   hu = flipud(h);
   zi = 1./z;
   for ka = 1:size(h,2)
      h(:,ka) = polyval(hu(:,ka),zi);
   end
   scales = [];
   p = 1:size(h,2); % use basis indices; there are no "poles" here
   L = 1;
else   
   p = uniform_over_ring_sector(opt.r1,opt.r2,opt.NumAtoms,opt.phi1,opt.phi2,C);
   if C
      pc = p(end);
      p = p(1:end-1);
   end
   p = p(:).'; % row vector
   Ir = imag(p)==0;
   nr0 = sum(Ir);
   [scales,L] = getScales(p,NumSamples); % nr+2*nc elements
   % duplicate scales for negative atoms.
   scales = [scales(1:nr0),scales(1:nr0),scales(nr0+1:end),scales(nr0+1:end)];
   
   h1 = 1./(z-p(Ir));
   h1a = 1./(z+p(Ir));
   Dr = z.^2-2*real(p(~Ir)).*z+abs(p(~Ir)).^2;
   Dra = z.^2+2*real(p(~Ir)).*z+abs(p(~Ir)).^2;
   h2 = 2*imag(p(~Ir))./Dr;
   h2a = 2*imag(p(~Ir))./Dra;
   h3 = 2*(z-real(p(~Ir)))./Dr;
   h3a = 2*(z+real(p(~Ir)))./Dra;
   h = [h1,h1a,h3,h2,h3a,h2a];
   h = h.*scales;
   p = [p(Ir),-p(Ir),p(~Ir),conj(p(~Ir)),-p(~Ir),conj(-p(~Ir))];
end

if C
   p = [p, 1];
   h = [h, ones(size(h,1),1)/L];
end

% delays
hh = z.^-(1:10)./(1:10);
h = [h,hh];
p = [p, zeros(1,10)];
