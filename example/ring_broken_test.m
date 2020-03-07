%% Initial Setup

load('Ring_Anorm_broken.mat');

%Optimum of previous reweighting iteration
BM =  opt.warm_start.bash_manager;
y  =  opt.warm_start.y;
x  =  BM.get_x(y);
opt.DG_tol = 1e-4;

grad_w = BM.gradient_full(y);
e_w = BM.get_error(y);

%Move onto optimum of new skewed face (preliminary bash, done in BB already)
[BM, y_bash] = BM.bash(y);
grad_bash  = BM.gradient_full(y_bash);
e_bash = BM.get_error(y_bash);
x_bash = BM.get_x(y_bash);

%check norm consistencies
[~, x_norm] = LMO_1d(x, opt.norm_type, opt.w);
[~, x_norm_bash] = LMO_1d(x_bash, opt.norm_type, opt.w);

opt.visualize_end = 1;
opt.visualize = 0;

%start from the previous point?
WARM_START = 1;
%Use weights=1 instead? (only when  WARM_START = 0)
SAME_WEIGHTS = 0;

%% Run BB algorithm
if ~WARM_START
    opt = rmfield(opt,'warm_start');
    
    if SAME_WEIGHTS && ~WARM_START
        opt.w.weights = ones(size(opt.w.weights));
    end
end



[x_final, S_final, c_final, run_log] = BB_operator(A, At, b, opt);
e_final = 0.5*norm(A(x_final)-b,2)^2 + 0.5*opt.delta*norm(x_final,2)^2;
 
%Notice the discrepency in norm, it does not equal 200
[~, x_norm_final] = LMO_1d(x_final, opt.norm_type,  opt.w);