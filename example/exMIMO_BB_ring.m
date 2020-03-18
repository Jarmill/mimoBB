load('ring_sys.mat')

SOLVE = 1;
FORWARD = 1;
DRAW = 1;

rng(77, 'twister')

ny = size(sys, 1);
nu = size(sys, 2);

Nf =  128;
%W = 0.1*ones(Nf, nu, ny);
W = 1e-2*ones(ny, nu, Nf);
%W = 1e-1*ones(ny, nu, Nf);
%W = [];

opt.FreqWeight = W;
%opt.FreqWeight = ones(Ns, ny);
%opt.Compare = 1;
opt.Compare = 0;
opt.tau = 200;
opt.RandomRounds = 10;
opt.ReweightRounds = 10;
opt.NumAtoms = 1e3;
opt.NormType = Inf;
opt.FormSystem = 1;
opt.FCFW = 1;
%opt.delta = 1;
if SOLVE
[out, out_random] = exTrace_BB(sys,Ns,SNR,bw,opt); % target cost: 83202.8

opt.FCFW = 0;
opt.RandomRounds = 0;
opt.ReweightRounds = 0;

[out_fw, out_random] = exTrace_BB(sys,Ns,SNR,bw,opt); % target cost: 83202.8
if opt.Compare
    utGenAnalysisPlots(out,sys) % quality analysis
end
end

if DRAW
figure(2)
FS = 18;
clf
% cm_viridis=viridis(m);
% colormap(cm_viridis)
%iopzmap(out.sys_out)
subplot(1,2,1)
hold on
th = linspace(0, 2*pi,  201);
scatter(real(out_random.poles_active), imag(out_random.poles_active), 200,'x')
plot(cos(th), sin(th), 'k') 
text(-0.4, 0, sprintf('Order:\nCost:'), 'Fontsize', FS)
text(0.1, 0, sprintf('%i\n%0.2e', out_random.system_order, out_random.cost), 'Fontsize', FS)
hold off
axis square
box off
title('Poles before reweighting', 'Fontsize', FS) 

xlabel('Re(z)')
ylabel('Im(z)')
xticks([-1,-0.5,0,0.5,1])
yticks([-1,-0.5,0,0.5,1])


subplot(1,2,2)
hold on
th = linspace(0, 2*pi,  201);
scatter(real(out.poles_active), imag(out.poles_active), 200, 'x')
plot(cos(th), sin(th), 'k') 
text(-0.4, 0, sprintf('Order:\nCost:'), 'Fontsize', FS)
text(0.1, 0, sprintf('%i\n%0.2e', out.system_order, out.cost), 'Fontsize', FS)
hold off
axis square
title('Poles after reweighting', 'Fontsize', FS) 
xlabel('Re(z)')
ylabel('Im(z)')
xticks([-1,-0.5,0,0.5,1])
yticks([-1,-0.5,0,0.5,1])

end
