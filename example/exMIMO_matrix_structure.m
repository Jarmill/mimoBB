Ns = 11; nu = 2; ny = 3; nx = 4;
SNR = 60; % signal_var/noise_var
rho = 0.8; % spectral radius

opt = sisoAtomOptions;
opt.ShowProgressPlot = true;
opt.IncludeConstant = true;
opt.r1 = rho; % is this cheating or reasonable prior knowledge

rng default
rng(1)
sys = utGenExampleSystem(rho,ny,nu,nx);
opt.tau = 62;
%opt.NumAtoms = 50;
opt.NumAtoms = 5;

s = sqrt(db2mag(SNR));
u = randn(Ns,nu);
y = lsim(sys,u); % note: zero IC
n = randn(Ns,ny);
yn = y;
for ky = 1:ny
   yn(:,ky) = y(:,ky) + norm(y(:,ky),2)/s*n(:,ky);
end
z = iddata(y,u,1);
zn = iddata(yn,u,1);

TargetCost = norm(yn-y)^2/2;

%% Toeplitz Matrix
% Tu = cell(1,nu);
% for ku = 1:nu
%    Tu{ku} = toeplitz([u(1,ku);zeros(Ns-1,1)],u(:,ku))';
% end
Tu = cell(1,nu);
F = cell(1,nu);
for ku = 1:nu    
   Tu{ku} = toeplitz(u(:,ku) ,[u(1,ku);zeros(Ns-1,1)]);
   F{ku}  = toeplitzmultaux(u(:,ku),  [u(1,ku);zeros(Ns-1,1)]);
   %matrix-vector toeplitz multiplication:
   %y2 = toeplitzmult2(F, x);
end

%poles
[ha,p] = createAtoms(Ns,opt);
np = size(ha, 2);
%ha: matrix of pole responses
%p:  poles
Tuu = cell2mat(Tu);
A = kron(eye(ny), Tuu);
kp = kron(eye(nu), ha);

Ap = kron(eye(ny), Tuu*kp);

Apt = kron(eye(ny ), kp'*Tuu');

%% Test the multiplications
%c = randn(nu, ny, np);
c = cell(nu, ny);
for i = 1:ny
    for j = 1:nu
        c{j,i} = randn(np, 1);
    end
end
%cr = reshape(permute(cell2mat(c), [2 1 3]), [], 1);
cr = reshape(cell2mat(c), [], 1);

%coefficients indexed by [y, u, p]
cc = permute(reshape(cr, np, nu, ny ), [3,2,1]);

Acref = Ap*cr;



%Test A*c
Ac = cell(1, ny);
for i = 1:ny
    Ac{i} = zeros( Ns, 1);
end
%Ac = zeros(Ns, ny);
c0 = [];
rc0 = [];
for i = 1:ny    
    Ac_curr = zeros(1, Ns);
    for j = 1:nu
        c_curr = c{j,i};
        rc_curr = ha*c_curr;
        urc_curr = toeplitzmult2(F{j}, rc_curr);
        Ac{i} = Ac{i} + urc_curr;
        c0 = [c0; c_curr];
        rc0 = [rc0; rc_curr];
    end    
end
Acr = reshape(cell2mat(Ac), [], 1);

Acr2 = mimo_A(cr, np, nu, ny, Ns, F, ha);

norm(Acr - Acref);
norm(Acr - Acr2);
%Acr = reshape(Ac, [], 1);

%Test A'b
Atbref = Apt*Acr;
b = Ac;

Atb = cell(nu, ny);
for i = 1:ny
    for j = 1:nu
        Atb{j,i} = zeros(np, 1);
    end
end

for j = 1:nu
    for i = 1:ny
        b_curr  = b{i};
        ub_curr = toeplitzmult2(conj(F{j}), b_curr);
        rub_curr = ha'*ub_curr;
        Atb{j,i} = rub_curr;
    end
end

Atbr = reshape(cell2mat(Atb), [], 1);

Atbr2 = mimo_At(Acr,np, nu, ny, Ns, F, ha);

norm(Atbr - Atbref)
norm(Atbr2 - Atbref)

% 
% Ac_diff = norm(Acr - Acref)
%c = randn(nu*ny*size(ha, 2), 1);

