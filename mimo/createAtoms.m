function [h,p_out,scales, groups, f] = createAtoms(NumSamples,opt)
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
FreqSample = opt.FreqSample;
p = uniform_over_ring_sector(opt.r1,opt.r2,opt.NumAtoms,opt.phi1,opt.phi2,C);

%generate the poles
[h,p_out, f, scales, groups] = createAtomsDict(p,FreqSample, C, NumSamples);

return

% delays
%hh = [diag(ones(1,10),-1);zeros(NumSamples-11,11)];
%hh = hh(:,1:end-1)./(1:10);
%h = [h,hh];
%p = [p, zeros(1,10)];


