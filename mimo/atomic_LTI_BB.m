%% MIMO IDENTIFICATION: L2, TIME DOMAIN
% Minimize prediction error trace L2 norm for MIMO LTI identification.
% We consider a partial fraction structure G(z) = \sum_i R_i/(z-p_i), where
% R is residue matrix of size [ny, nu]. For mechanical structures, R is a
% dense matrix. But in general, we impose low cardinality of vec(R).
%
%Fully Corrective Bag-and-bash implementation
%
%Use sum-of-norms regularization for MIMO system


% Input fields:
% u     = Input Toeplitz matrix
% ym     = Measured signal (Output)
% tau    = Bound on atomic norm of the solution
% t_max  = Maximum iteration number
% p      = Random poles
% m_pole = Number of poles to try @ each iteration
% h0     = (Optional) Initial value of the solution

%% Options
% 1a. MLE: Trace vs Det
% 1b. Worst case: Inf-norm
% 2. IC
% 3. Missing data
% 4. time and freq data
% 5. Constraints: Damping, Natural frequency, Max Overshoot, Bandwidth, DC
%    gain, frequency weight, frequency bounds

function out = atomic_LTI_BB(In,opt)
y = In.ym;
u = In.u;


tau = opt.tau;
t_max = opt.MaxIter;
[Ns,ny] = size(y);
nu = size(u,2);


%Toeplitz the input for fast multiplication
F = cell(1,nu);
for ku = 1:nu    
   %Tu{ku} = toeplitz(u(:,ku) ,[u(1,ku);zeros(Ns-1,1)]);
   F{ku}  = toeplitzmultaux(u(:,ku),  [u(1,ku);zeros(Ns-1,1)]);
   %matrix-vector toeplitz multiplication:
   %y2 = toeplitzmult2(F, u);
end

%Process the poles into groups for sum-of-norms regularization
p =  In.PoleArray;
ha = In.ImpRespArray;
np = size(p, 2);

g = In.PoleGroups;
g_hot = ind2vec(g)'; %one-hot encoding of groups
g_offset = (0: (ny*nu-1))*np;
Ngroups = max(g);

w.groups = cell(Ngroups, 1);
w.weights = zeros(Ngroups, 1);
w.order = zeros(Ngroups, 1);
for gi = 1:Ngroups
    i_curr = find(g_hot(:, gi));
    %penalize the use of second-order poles
    if length(i_curr) == 1
        w.weights(gi) = 1;
        w.order(gi) = 1;
    else
        %w.weights(gi) = sqrt(2); %sqrt(2)?
        w.weights(gi) = 2; %sqrt(2)?
        w.order(gi) = 2;
    end
    
    ind_gi = i_curr+g_offset;
    w.groups{gi} = ind_gi(:);    
end

%Formulate the least squares operators and paramters
%system to output wrt. input and its adjoint
%A  = @(x) mimo_A(x, np, nu, ny, Ns, F, ha);
%At = @(r) mimo_At(r,np, nu, ny, Ns, F, ha);
A  = @(x) mimo_A2(x, np, nu, ny, Ns, F, ha);
At = @(r) mimo_At2(r,np, nu, ny, Ns, F, ha);

b = reshape(y, Ns*ny, 1);

BB_opt.num_var = nu*ny*np;
BB_opt.tau = tau;
BB_opt.w = w;
BB_opt.delta = 0;
%BB_opt.norm_type = 2;
BB_opt.norm_type = Inf;
BB_opt.is_complex = 0;

tic
%Run the optimization routine
[x_final, S_final, c_final, run_log] = BB_operator(A, At, b, BB_opt);

toc


out.Coeff0 = x_final;
out.Atoms = S_final;
out.AtomCoeff = c_final;
out.run_log = run_log;

Ax = A(x_final);
%h_all = A(x_final)

out.y = reshape(Ax, Ns, ny);

%determine output
x_coeff = reshape(x_final, np, nu, ny);

out.Coeff = cell(nu, ny);
out.h = cell(nu, ny);

for j = 1:nu
    for i = 1:ny
        x_curr = x_coeff(:, j, i);
        out.Coeff{j, i} = x_curr;
        out.h{j, i} = ha*x_curr;
    end
end

out.iter = length(run_log.time);
out.cost = norm(Ax - b)^2/2;

out.group_active =  zeros(Ngroups, 1);
for gi = 1:Ngroups
    g_curr = w.groups{gi};    
    
    if any(x_final(g_curr))
        out.group_active(gi) = w.order(gi);
    end
    
end

out.system_order = sum(out.group_active);

%h_svd = svd(hankel_mo(h(2:min(1000,N),:)'));
%out.svd = h_svd(1:min(15,length(h_svd)))';
%assignin('base','out',out);

%out.error = norm(In.h-out.h)/norm(In.h);
end