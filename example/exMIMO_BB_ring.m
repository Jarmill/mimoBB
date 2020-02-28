load('ring_sys.mat')

rng(77, 'twister')

ny = size(sys, 1);
nu = size(sys, 2);

Nf =  128;
%W = 0.1*ones(Nf, nu, ny);
W = 0.01*ones(ny, nu, Nf);

opt.FreqWeight = W;
%opt.FreqWeight = ones(Ns, ny);
%opt.Compare = 1;
opt.Compare = 0;
opt.tau =  700;
%opts.ReweightRounds = 10;
out = exTrace_BB(sys,Ns,SNR,bw,opt); % target cost: 83202.8

if opt.Compare
    utGenAnalysisPlots(out,sys) % quality analysis
end