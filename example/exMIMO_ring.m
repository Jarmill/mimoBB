% outputs.

addpath ..
close all

load('ring_sys.mat')

Ns = 1000; nu = 3; ny = 2; nx = 4;
SNR = 60; % signal_var/noise_var
opt = sisoAtomOptions;
opt.ShowProgressPlot = true;
opt.IncludeConstant = true;
opt.ModelType = "ss"; 
opt.SearchMethod = "grad";
opt.Type = "rational"; 
opt.IsDet = false;
if opt.Type=="rational"
   opt.MaxIterTrace = 400;
   opt.NumAtoms = 500;
else   
   % not implemented yet
   %
end

rho = 0.8; % spectral radius
opt.r1 = rho; % is this cheating or reasonable prior knowledge

%% Example 1:
rng default
rng(1)
%sys = utGenExampleSystem(rho,ny,nu,nx);
%opt.tau = 62;
opt.tau = 200;
opt.NumAtoms = 200;
opt.MaxIterTrace = 800;
opt.ShowProgressPlot = false;
opt.SearchMethod = "grad";
opt.Randomize = true;
opt.AwayStepOnNeedBasis = true;
out = exTrace(sys,Ns,SNR,opt); % target cost: 83202.8
utGenAnalysisPlots(out,sys) % quality analysis
