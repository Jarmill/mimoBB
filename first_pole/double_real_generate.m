%generate a system with two real poles

% Ts = 0.05;
% Ts = 1;
Ts = 0.1;
k = 6;
G = zpk([0, 0.9], [0.8, -0.5], -k, Ts);

%generate the step-response data
% Thorizon = 25;
Thorizon = 10;
t = (0:Ts:Thorizon)';

%impulse response
% y=impulse(G,t);
% u=zeros(size(y));
% u(1) = 1;
% y=step(G,t);
% t = 0:0.04:8;  % 201 points

% prng system
u = 2*rand(size(t))-1;
y = lsim(G, u, t);


opt = mimoAtomOptions;
opt.SampleTime = Ts;
% opt.phi1 = 0;
% opt.phi2 = 0;
opt.tau = 5;

[out, out_random] = atomic_LTI_iteration(y, u, opt);

% figure(1)
% clf
% plot(t, y)
plot(t, [y, out.y])