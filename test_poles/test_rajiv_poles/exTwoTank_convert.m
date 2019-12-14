% Two Tank Data System Identification
% Model format: PVNARX with NL = []
rng default
rng(3)
load twotankdata; % ships with System Identification Toolbox
Ts = 0.2;
z = iddata(y,u,Ts);
z = z(1:2:end);
ze = z(1:500); ze0 = ze;
zv = z(576:end);
Ns = size(ze,1);
ze.Ts = 1;
ze.Tstart = 0;
zv.Ts = 1;
zv.Tstart = 0; % for simulink

% narx orders
na = 2; nb = 2; nk = 1;
% create template model (linear in regressor form)
m0 = idnlarx([na nb nk],[]); % [] means f(.) is linear
%Rp = polyreg(m0,'MaxPower',2,'CrossTerm','off');
%m0 = addreg(m0, Rp);
m1 = nlarx(ze,m0); % benchmark
%m2 = nlarx(ze0,m0);
R = getreg(m0,'all',ze); 
R = [R,ones(Ns,1)]; 

C = polyreg(m0,'max',2,'cross','on');
m2 = nlarx(ze0,[2 2 1],'sigmoid','custom',C);

% Unlike IDNLARX model, we don't have a structural element to account for
% offsets. Hence to account for offsets, include a column vector of ones in
% R. For now, assume offset is zero. For theta-dynamics, always include a
% real pole at 1 to account for constant portion:
% Theta(i) = [A0 + \sum a_i*Ai], where A0 is a constant scalar.

ye0 = ze.y; ue = ze.u;
ye = [0; ye0(1:end-1)];
Tu = toeplitz([ue(1);zeros(Ns-1,1)],ue)';
Ty = toeplitz([ye(1);zeros(Ns-1,1)],ye)';
%Tuy = toeplitz([ue(1).*ye(1);zeros(Ns-1,1)],ue.*ye)';
%In = struct('T',{{Tu,Ty}},'ym',ye,'tau',struct('tauAtom',0.0445),...
%   't_max',50,'k',100,'idx_avail',1:Ns,'Regmat',R,'alpha',0.6);
%In.h0 = out.h;

%% put the system in standard form for BB identification


%diagonal matrix of regressors
num_regs = size(R, 2);
Rd = reshape(R, [], 1);
%Rdd = sparse(kron(ones(1, num_regs), 1:length(R)), 1:length(Rd), Rd);
Rdd = spdiags(Rd, 0, length(Rd), length(Rd));

%Toeplitz matrix blocks
Tmt = sparse([Tu Ty]);
num_reg = size(R, 2);
Tmt_kron = kron(eye(num_reg), Tmt);
%Tmt_kron = kron(ones(1, num_reg), Tmt);

%standard atomic standard form
es = kron(ones(5, 1), speye(size(R, 1)));

%Hopefully this is the right linear transformation for the system
%after validating sysid, then work on taking advantage of system structure
%for efficient gradient/stepsize evaluation. Use kronecker forms and
%toeplitz fft multiplications.
A = es' * Rdd * Tmt_kron;


b = ye;
tau = 0.0445;
delta = 1e-5;

%system description
norm_type = 'poles';
w.n = 1000;
w.rho1 = 0.5;
w.rho2 = 0.8;
w.max_angle = pi/3;

%set up the different systems
sys_index = reshape(1:size(Tmt_kron, 2), [], 2*num_reg);
w.systems = num2cell(sys_index, 1);

[x_final, S_final, c_final, run_log] = BB_1d(A, b, tau, delta, norm_type, w);