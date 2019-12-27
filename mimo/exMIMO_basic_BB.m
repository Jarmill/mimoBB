% 2=outut, 3-input MIMO state-space (shared pole) LTI identification under
% atomic norm constraints. No extra demans are placed on the model other
% than fitting the measured output data in the trace sense, that is minimize
% trace(EE') under atomic norm constraint on each impulse response. E is
% Ns-by-ny matrix, where Ns := number of observations, ny := number of
% outputs.

addpath ..
close all

Ns = 1000; nu = 3; ny = 2; nx = 4;
SNR = 60; % signal_var/noise_var
opt = sisoAtomOptions;
opt.ShowProgressPlot = true;
opt.IncludeConstant = true;
opt.ModelType = "ss"; 
opt.SearchMethod = "grad";
opt.Type = "rational"; 
if opt.Type=="rational"
   opt.MaxIter = 400;
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
sys = utGenExampleSystem(rho,ny,nu,nx);


%opt.tau = 62;
opt.tau = 400;


opt.NumAtoms = 200;
opt.MaxIter = 800;
opt.ShowProgressPlot = false;
opt.SearchMethod = "grad";
opt.Randomize = true;
opt.AwayStepOnNeedBasis = true;
out = exTrace_BB(sys,Ns,SNR,opt); % target cost: 83202.8
utGenAnalysisPlots(out,sys) % quality analysis
% 
% 
% %% Example 2: repeated real pole (hard to get good results)
% rng default
% rng(2)
% sys = utGenExampleSystem(rho,ny,nu,nx);
% opt.tau = 32;
% opt.NumAtoms = 300;
% opt.MaxIter = 800;
% opt.ShowProgressPlot = false;
% opt.SearchMethod = "grad";
% opt.Randomize = false;
% opt.AwayStepOnNeedBasis = true;
% out = exTrace(sys,Ns,SNR,opt); 
% utGenAnalysisPlots(out,sys) 
% 
% %% Example 3: somewhat better results here
% rng default
% rng(3)
% sys = utGenExampleSystem(rho,ny,nu,nx);
% opt.tau = 12;
% opt.NumAtoms = 800;
% opt.MaxIter = 800;
% opt.ShowProgressPlot = false;
% opt.SearchMethod = "grad";
% opt.Randomize = false;
% opt.AwayStepOnNeedBasis = true;
% out = exTrace(sys,Ns,SNR,opt);  % 9925.85
% utGenAnalysisPlots(out,sys) 
% 
