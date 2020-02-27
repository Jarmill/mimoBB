%% SISO Challenge
% 1-output, 1-input SISO state-space (shared pole) LTI identification under
% atomic norm constraints. No extra demands are placed on the model other
% than fitting the measured output data in the trace sense, that is minimize
% trace(EE') under atomic norm constraint on impulse response. E is
% Ns-by-1 vector, where Ns := number of observations.

addpath ..
close all

Ns = 500; nu = 1; ny = 1; nx = 5;
SNR = 0; % signal_var/noise_var, in dB
opt = sisoAtomOptions;
opt.IncludeConstant = false;
opt.ModelType = "ss"; 
opt.Alpha = 0.75;
opt.SearchMethod = "grad";
opt.AdaptiveTau = false;

rho = 0.8; % spectral radius
opt.r1 = rho; % is this cheating or reasonable prior knowledge

rng default
rng(3)
sys = utGenExampleSystem(rho,ny,nu,nx);
opt.tau = 8;
opt.Type = "rational";
opt.NumAtoms = 1000;

opt.FW = "FA";
opt.phi2 = pi*0.1;
opt.ShowProgressPlot = false;
opt.SearchMethod = "grad";
opt.Randomize = false;
opt.AwayStepOnNeedBasis = true;
opt.IsDet = false;

if opt.IsDet
   opt.MaxIterTrace = 100;
   opt.MaxIterDet = 15;
else
   opt.MaxIterTrace = 500;
end
out = exDet(sys,Ns,SNR,0.1,opt); 
%utGenAnalysisPlots(out,sys) 

figure; 
k = 1;
L = floor(Ns/2);
H=hankel(out.h(2:L,k),out.h(L:end,k));bar(svd(H));xlim([0 10]);

figure
c=exp(1i*2*pi*(0:2000)/2000);
[ha,p] = createAtoms(Ns,opt);
plot(real(p),imag(p),'.',real(c),imag(c),'k')

% true NV
n=out.Results.zn(1:1000).y-out.Results.z(1:1000).y;NV0=n'*n/1000;

figure
pzmap(out.Results.sys,out.Results.sys_ssest_td,out.Results.sysr)
