function [x_final, run_log] = BB_forward(data, b, opt)
%simple implementation of forward frank wolfe
%same arguments as BB_operator

tau = opt.tau;
delta = opt.delta;
w = opt.w;
N = opt.num_var;
norm_type = opt.norm_type;
DG_tol = opt.DG_tol;

if isnumeric(data)
    %can probably do this more elegantly
    A0 = data;
    At = @(r) A0'*r;
    A  = @(x) A0*x;
    A_matrix = 1;
else
    At = data.At;
    A  = data.A;
    A_matrix = 0;
end


terminate = 0;
k = 1;

res = Ax - b;
grad = At(res) + delta*x;
[n, nlist] = anorm_1d(x,  opt.norm_type, opt.w);

error_orig = 0.5*norm(Ax-b)^2 + 0.5*delta*norm(x)^2;
error_old = error_orig;
error_gap = Inf;


run_log.error_list  = [];
run_log.atomic_norm = [];
%run_log.atomic_norm_true = [];
%run_log.cardinality = [];
%run_log.num_attempted = [];
%run_log.num_survived = [];
run_log.duality_gap = [];
run_log.time = [];
end

tic
while ~terminate  
    %[BM, S_bag, DG] = BM.bag_atoms(grad, x, N_bag);
    [a, n] = LMO_1d(-grad, norm_type, w);
    wdir = (tau*a - x);
    DG = -grad'*wdir;
    %if the bag is empty, then no more atoms can be added to the system
    %this is a generalized termination condition
    %if isempty(S_bag) || abs(error_gap) < 1e-7
    % k > 1 && (rank(full([S_new S_bag])) ~= (size(S_new, 2)  + 1))
    if DG > DG_tol
        Awdir = A(wdir);
        
        alpha_top = -DG;
        alpha_bottom = norm(Awdir,2)^2 + delta*norm(wdir,2)^2;
        alpha = max(0, min(alpha_max, alpha));        
        
        x_new = x + alpha * wdir;
        Ax_new = A(x_new);
        grad_new = At(Ax_new - b);
        error_new = norm(Ax_new, 2)^2 + delta*norm(x_new, 2)^2;
        n_new = anorm_1d(x_new, norm_type, w);
    else
        terminate = 1;
        x_new = x; 
        error_new = error;
        n_new = n;
    end
    
    error_gap = error - error_new;
    
    run_log.error_list(k) = error_new;
    run_log.atomic_norm(k) = n_new;
    run_log.duality_gap(k) = DG;
    run_log.time(k) = toc;
    
    x = x_new;
    Ax = Ax_new;
    grad = grad_new;
    n = n_new;
    
end

toc