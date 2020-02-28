% 2-output, 3-input MIMO state-space (shared pole) LTI identification under
% atomic norm constraints. No extra demands are placed on the model other
% than fitting the measured output data in the trace sense, that is minimize
% trace(EE') under atomic norm constraint on each impulse response. E is
% Ns-by-ny matrix, where Ns := number of observations, ny := number of
% outputs.

addpath ..
close all
warning off Ident:estimation:invalidFocusOption2

Ns = 200; nu = 3; ny = 2; nx = 5;
SNR = 0; % signal_var/noise_var
opt = sisoAtomOptions;
opt.ShowProgressPlot = true;
opt.IncludeConstant = true;
opt.ModelType = "ss"; 
opt.SearchMethod = "grad";
opt.Type = "rational"; 

%opt.Freq = "time";

if opt.Type=="rational"
   opt.MaxIterTrace = 400;
   opt.NumAtoms = 500;
else   
   % not implemented yet
   %
end

rho = 0.97; % spectral radius
opt.r1 = rho; % is this cheating or reasonable prior knowledge
bw = 0.1;

%% Example 1:
rng default
rng(1)
%{
sys = utGenExampleSystem(rho,ny,nu,nx);
[U,S,V]=svd(sys.A); 
S(1,1)=0.999;
%S(1,1)=0.999; S(2,2)=0.985;  S(3,3)=0.975; 
sys.A=U*S*V';
%}
load testSys1

%opt.tau = 62;
%opt.tau = 400;
%opt.tau = 800;
opt.tau = 4000;


opt.NumAtoms = 1000;
opt.MaxIterTrace = 800;
opt.ShowProgressPlot = false;
opt.SearchMethod = "grad";
opt.Randomize = true;
opt.AwayStepOnNeedBasis = true;
opt.ReweightRounds = 10;


%W = ones(Ns, nu, ny);
W = ones(ny, nu, Ns);
%opt.FreqWeight = W;
opt.FreqWeight = [];

out = exTrace_BB(sys,Ns,SNR,bw,opt); % target cost: 83202.8
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
