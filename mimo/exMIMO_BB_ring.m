load('ring_sys.mat')

rng(77, 'twister')

ny = size(sys, 1);
nu = size(sys, 2);
%W = ones(Ns, nu, ny);
%opt.FreqWeight = W;
%opt.FreqWeight = ones(Ns, ny);
opt.Compare = 1;
%opt.tau =  600;
%opts.ReweightRounds = 10;
out = exTrace_BB(sys,Ns,SNR,bw,opt); % target cost: 83202.8
utGenAnalysisPlots(out,sys) % quality analysis