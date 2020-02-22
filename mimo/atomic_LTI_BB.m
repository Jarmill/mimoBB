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




function out = atomic_LTI_BB(In,opt)
y = In.ym;
u = In.u;

tau = opt.tau;
t_max = opt.MaxIterTrace;
[Ns,ny] = size(y);
nu = size(u,2);

Y = fft(y, Ns, 1);
U = fft(u, Ns, 1);

%Toeplitz the input for fast multiplication
F = cell(1,nu);
for ku = 1:nu    
   %Tu{ku} = toeplitz(u(:,ku) ,[u(1,ku);zeros(Ns-1,1)]);
   F{ku}  = toeplitzmultaux(u(:,ku),  [u(1,ku);zeros(Ns-1,1)]);
   %matrix-vector toeplitz multiplication:
   %y2 = toeplitzmult2(F, u);
end



%Process the poles into groups for sum-of-norms regularization
p  = In.PoleArray;
ha = In.ImpRespArray;
f  = In.FreqRespArray;
W  = In.FreqWeight;
np = size(p, 2);

g = In.PoleGroups;
g_hot = ind2vec(g)'; %one-hot encoding of groups
%g_offset = (0: (ny*nu-1))*np;
%g_offset = 0:(nu*ny-1);
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
    
    ind_gi = (i_curr-1)*ny*nu+(1:ny*nu);
    w.groups{gi} = ind_gi(:);    
end

%Formulate the least squares operators and paramters
%system to output wrt. input and its adjoint

b_time = reshape(y, Ns*ny, 1);

%Frequency response penalization
Wdim = length(size(W));


%check the b_freq calculations

if Wdim == 3
    %IO (weighting function for each input/output pair)
    A  = @(x) mimo_io_A(x, np, nu, ny, Ns, F, ha, f, U, W);
    At = @(r) mimo_io_At(r,np, nu, ny, Ns, F, ha, f, U, W);
    b_freq = complex_unfold(kron(reshape(permute(Y, [2,1]), [], 1), ones(nu, 1)), 1);
elseif Wdim == 2
    %Output (weighting function for each output)
    A  = @(x) mimo_output_A(x, np, nu, ny, Ns, F, ha, f, U, W);
    At = @(r) mimo_output_At(r,np, nu, ny, Ns, F, ha, f, U, W);
    b_freq = complex_unfold(reshape(permute(Y, [2,1]), [], 1));
else   
    %Time (no frequency penalization)
    A  = @(x) mimo_A2(x, np, nu, ny, Ns, F, ha);
    At = @(r) mimo_At2(r,np, nu, ny, Ns, F, ha);
    b_freq = [];
end




%b_freq = squeeze(reshape(G_ref, [], 1, 1));

%b_freq = repmat(Y, 1, 1, nu);
%

%I really hope this works
%want to copy Y, slicing by indices
%probably easier to write it out or keep everything as arrays
%or just let reshape take care of everything
%b_freq = kron(reshape(permute(Y, [2,1]), [], 1), [1;1;1]);
%b_freq = squeeze(reshape(repmat(permute(Y, [2,1]), 2,1,1), [], 1));
%
%b_freq = squeeze(reshape(b_freq, 1, 1 ,[] ));


b = [b_time; b_freq];

%Weighting function?



%something about reshaping?

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


%output from data
Ax = A(x_final);
%h_all = A(x_final)

%time
y_time = b(1:(Ns*ny));
out.y = reshape(y_time, Ns, ny);

%Frequency
y_freq_real = b(Ns*ny + 1:end);
y_freq = complex_fold(y_freq_real, 1);

if Wdim == 3
    out.f = reshape(y_freq, Ns, nu, ny);
elseif Wdim == 2
    out.f = reshape(y_freq, Ns, ny);
else
    out.f = [];
end
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


%reweighting?
delta = 1e-4;
out.group_active =  zeros(Ngroups, 1);
out.weights_new = [];
for gi = 1:Ngroups
    g_curr = w.groups{gi};    
    
    if any(x_final(g_curr))
        out.group_active(end+1) = gi;
        out.weights_new(end+1) = 1/(delta + max(x_final(g_curr)));
        group_poles(gi) = w.order(gi);
    end
    
end


out.system_order = sum(out.group_active);

gaa = any(In.PoleGroups == out.group_active, 1);
out.poles_active = In.PoleArray(gaa)';


%h_svd = svd(hankel_mo(h(2:min(1000,N),:)'));
%out.svd = h_svd(1:min(15,length(h_svd)))';
%assignin('base','out',out);

%out.error = norm(In.h-out.h)/norm(In.h);
end