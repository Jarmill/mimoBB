load('ring_sys.mat')

rng(77, 'twister')

ny = size(sys, 1);
nu = size(sys, 2);

Nf =  128;
%W = 0.1*ones(Nf, nu, ny);
W = 0.05*ones(ny, nu, Nf);

opt.FreqWeight = W;
%opt.FreqWeight = ones(Ns, ny);
%opt.Compare = 1;
opt.Compare = 0;
opt.tau =  600;
%opt.ReweightRounds = 6;
opt.RandomRounds = 15;
out = exTrace_BB(sys,Ns,SNR,bw,opt); % target cost: 83202.8

figure(2)
clf
iopzmap(out.sys_out)

if opt.Compare
    utGenAnalysisPlots(out,sys) % quality analysis
end