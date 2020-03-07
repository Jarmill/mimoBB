load('ring_sys.mat')

rng(77, 'twister')

ny = size(sys, 1);
nu = size(sys, 2);

Nf =  128;
%W = 0.1*ones(Nf, nu, ny);
%W = 0.05*ones(ny, nu, Nf);
W = 1e-1*ones(ny, nu, Nf);
%W = [];

opt.FreqWeight = W;
%opt.FreqWeight = ones(Ns, ny);
%opt.Compare = 1;
opt.Compare = 0;
opt.tau = 200;
opt.RandomRounds = 20;
opt.ReweightRounds = 25;
opt.NumAtoms = 1e3;
opt.NormType = Inf;
out = exTrace_BB(sys,Ns,SNR,bw,opt); % target cost: 83202.8

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

if opt.Compare
    utGenAnalysisPlots(out,sys) % quality analysis
end