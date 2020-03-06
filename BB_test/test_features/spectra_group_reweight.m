load 'spectra'

tau_max = 100;
opt = struct;
opt.num_var = size(NIR, 2);
opt.tau = tau_max;
opt.delta = 0;
opt.norm_type = 2;
opt.is_complex = 0;
opt.export_warm_start = 1;

w0 = [1,1,1,1,1];

opt.w = struct;
opt.w.groups = {[1:80], [81:160], [161:240], [241:320], [321:401]};
opt.w.weights = w0;
Ng = length(opt.w.weights);

%BB_1d
Af  = @(x) data_A(NIR, x);
Atf = @(x) data_At(NIR, x);

[x_g1, S_g1, c_g1, run_g1] = BB_operator(Af, Atf, octane, opt);

epsilon =  1e-4;

anorm =  0;

ReweightRounds = 10;
weights_list = cell(1+ReweightRounds, 1);
weights_list{1} = opt.w.weights;
vnorm_list = cell(ReweightRounds, 1);
error_list = zeros(ReweightRounds, 1);
x_list = zeros(opt.num_var, ReweightRounds);
for k = 1:ReweightRounds
    [x_g1, S_g1, c_g1, run_g1] = BB_operator(Af, Atf, octane, opt);
    x_list(:,  k) = x_g1;
    error_list(k) = run_g1.error_list(end);
    weights_new = [];
    for i = 1:Ng
        vnorm(i)= norm(x_g1(opt.w.groups{i}), opt.norm_type);
        anorm = anorm + opt.w.weights(i)* vnorm(i);

        weights_new(i) = w0(i)/(epsilon + vnorm(i));
    end
    
    weights_new_norm = weights_new * opt.tau  / (vnorm*weights_new');
    weights_list{k+1} = weights_new_norm;
    vnorm_list{k} = vnorm;
    opt.w.weights = weights_new_norm;
end


% %reweight 1
% x_ind = find(x_g1);
% x_active = x_g1(x_ind);
% w0 = 1./(abs(x_active) + epsilon);
% w = w0 *  (tau_max/ (w0'*x_active));
% 
% NIR_active  = NIR(:, x_ind);
% 
% Ar  = @(x) data_A(NIR_active, x);
% Atr = @(x) data_At(NIR_active, x);
% 
% opt.w = w;
% opt.num_var = size(NIR_active, 2);
% [x_r, S_r, c_r, run_r] = BB_operator(Ar, Atr, octane, opt);
% 
% %reweight2
% x_ind1 = find(x_r);
% x_active1 = x_r(x_ind1);
% w_active1 = w(x_ind1);
% w01 = 1./(abs(x_active1) + epsilon);
% w1 = w01 *  (tau_max/ (w01'*abs(x_active1)));
% %w1 =  w(x_ind1);
% 
% NIR_active1  = NIR_active(:, x_ind1);
% 
% Ar1  = @(x) data_A(NIR_active1, x);
% Atr1 = @(x) data_At(NIR_active1, x);
% 
% opt.w = w1;
% opt.num_var = size(NIR_active1, 2);
% [x_r1, S_r1, c_r1, run_r1] = BB_operator(Ar1, Atr1, octane, opt);
% 
% %x_pw = ANEN_helpers(NIR, octane, tau, delta);
% 
% x_r_out = zeros(size(x_g1));
% x_r_out(x_ind) = x_r;
% 
% x_r_out1 = zeros(size(x_g1));
% x_r_out1(x_ind(x_ind1)) = x_r1;
% 
figure(25)
clf
hold on
stem(x_list, '.')
% stem(x_r_out, '.')
% stem(x_r_out1, '.')
% %stem(x_max, '.')
% %stem(x_p, '.')
% %stem(x_en, '.')
% %stem(x_delta, '.')
% %stem(x_pw, '.')
% %legend('Constrained min (affine)', 'L2 min')
hold off
% 
%error_g1 = run_g1.warm_start.bash_manager.get_error(run_g1.warm_start.y);
% error_r  = run_r.warm_start.bash_manager.get_error(run_r.warm_start.y);
% error_r1  = run_r1.warm_start.bash_manager.get_error(run_r1.warm_start.y);

function Ax = data_A(A, x)
    Ax = A*x;
end

function Atx = data_At(A, x)     
    Atx = A'*x;    
end
