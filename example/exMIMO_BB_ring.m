load('ring_sys.mat')

SOLVE = 1;
FORWARD = 1;
DRAW = 1;

rng(77, 'twister')

ny = size(sys, 1);
nu = size(sys, 2);

Nf =  128;
wc = 2.5;
k = 0.1;
% W = @(f) k ./((2*pi*1.0j.*f)/wc +1);
% W = 0;
%corner frequency
% wc = 2;


% W = 1e-2*ones(ny, nu, Nf);
W = zeros(ny, nu, Nf);


%W = 1e-1*ones(ny, nu, Nf);
%W = [];
opt = mimoAtomOptions;
opt.FreqWeight = W;
%opt.FreqWeight = ones(Ns, ny);
%opt.Compare = 1;
opt.Compare = 0;
% opt.tau = 200;
% opt.tau = 250;
% opt.tau = 300;
% opt.tau = 500;
opt.tau = 600;
% opt.tau = 2000;
opt.RandomRounds = 50;
opt.ReweightRounds = 20;
opt.NumAtoms = 4e2;
opt.NormType = Inf;
opt.FormSystem = 1;
opt.FCFW = 1;


%sector bounds
opt.r1 = 0.98;
opt.phi2 = 0.5;
%opt.delta = 1;
if SOLVE
[out, out_random] = exTrace_BB(sys,Ns,SNR,bw,opt); % target cost: 83202.8

% opt.FCFW = 0;
% % opt.RandomRounds = 0;
% opt.ReweightRounds = 0;

% [out_fw, out_random] = exTrace_BB(sys,Ns,SNR,bw,opt); % target cost: 83202.8
% if opt.Compare
%     utGenAnalysisPlots(out,sys) % quality analysis
% end

sys_all = 0;
for i = 1:length(out.sys_modes)
    sys_all = sys_all + out.sys_modes{i};
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

figure(3)
clf
subplot(2, 1, 1)
hold on
plot(out.y_orig(:, 1))
plot(out.y(:, 1))
xlabel('time')
ylabel('y1')


subplot(2, 1, 2)
hold on
plot(out.y_orig(:, 2))
plot(out.y(:, 2))
xlabel('time')
ylabel('y2')


end
