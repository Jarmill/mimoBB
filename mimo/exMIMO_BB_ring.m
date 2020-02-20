load('ring_sys.mat')

rng(77, 'twister')

out = exTrace_BB(sys,Ns,SNR,bw,opt); % target cost: 83202.8
%utGenAnalysisPlots(out,sys) % quality analysis