% 2-output, 3-input MIMO state-space (shared pole) LTI identification under
% atomic norm constraints. No extra demands are placed on the model other
% than fitting the measured output data in the trace sense, that is minimize
% trace(EE') under atomic norm constraint on each impulse response. E is
% Ns-by-ny matrix, where Ns := number of observations, ny := number of
% outputs.

addpath ..
close all

Ns = 1000; nu = 3; ny = 2; nx = 5;
%Ns = 5; nu = 1; ny = 1; nx = 5;
SNR = 30; % signal_var/noise_var
opt = mimoAtomOptions;
opt.ShowProgressPlot = true;
opt.IncludeConstant = true;
opt.ModelType = "ss"; 
opt.SearchMethod = "grad";
opt.Type = "rational"; 

opt.FreqWeight = [];

if opt.Type=="rational"
   opt.MaxIterTrace = 400;
   opt.NumAtoms = 600;
else   
   % not implemented yet
   %
end

%rho = 0.97; % spectral radius
rho = 0.6; % spectral radius
opt.r1 = rho; % is this cheating or reasonable prior knowledge
bw = 0.1;

%% Example 1:
rng default
rng(1)
sys = utGenExampleSystem(rho,ny,nu,nx);


%opt.tau = 62;
%opt.tau = 400;
%opt.tau = 350; %this is a stress test
opt.tau = 200;

%opt.tau = 100;


opt.NumAtoms = 800;
opt.MaxIterTrace = 800;
opt.ShowProgressPlot = false;
%opt.SearchMethod = "grad";
%opt.Randomize = true;


%W = ones(Ns, nu, ny);
opt.FreqWeight = [];
opt.RandomRounds = 10;
opt.ReweightRounds = 10;

opt.Compare = 0;
out = exTrace_BB(sys,Ns,SNR,bw,opt); % target cost: 83202.8
%utGenAnalysisPlots(out,sys) % quality analysis

figure(2)
clf
%iopzmap(out.sys_out)
hold on
th = linspace(0, 2*pi,  201);
scatter(real(out.poles_active), imag(out.poles_active), 200, 'x')
plot(cos(th), sin(th), 'k') 
hold off
axis square
title('Poles of Identified System') 
xlabel('Re(z)')
ylabel('Im(z)')


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
