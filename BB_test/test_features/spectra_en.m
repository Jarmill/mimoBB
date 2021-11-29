load 'spectra'

%tau_max = 4;
%tau_max = 3;
%tau_max = 2.5;
%tau_max = 2;
%tau_max = 1.9;
%tau_max = 1.7;
%tau_max = 1.5;
%tau_max = 1.4;
%tau_max = 1.37;
%tau_max = 1.3;
%tau_max = 0.8;
%tau = 100;
%delta = 1.5;

tau_max = 100;

opt.num_var = size(NIR, 2);
opt.tau = tau_max;
opt.delta = 0;
opt.norm_type = 1;
opt.is_complex = 0;
opt.w = 1;
opt.export_warm_start = 1;


%x_en = ANEN_helpers(NIR, octane, tau, delta);
%x_no_delta = NIR \ octane;
%x_delta = spline_tikh_solve(NIR, octane, delta);
%K = NIR' * NIR;
%Atb= NIR'*octane;
%x_ols = K \ Atb;

%x_en = ANEN_affine(NIR, octane, tau, delta);
%x_en = ANEN_affine_tcp(NIR, octane, tau, delta);
%[x_l1, S_l1, c_l1] = BB_EN(NIR, octane, tau, delta, 1);
%[x_max, S_max, c_max] = BB_L_INF(NIR, octane, tau_max, delta, 1);
%[x_p, S_p, c_p] = BB_L_P(NIR, octane, tau_p, delta, p, 1);

%BB_1d
[x_l1, S_l1, c_l1, run_l1] = BB_operator(NIR, octane, opt);
%[x_max, S_max, c_max] = BB_1d(NIR, octane, tau_max, delta, Inf);
%[x_p, S_p, c_p] = BB_1d(NIR, octane, tau_p, delta, p);

epsilon =  1e-4;

%reweight 1
x_ind = find(x_l1);
x_active = x_l1(x_ind);
w0 = 1./(abs(x_active) + epsilon);
w = w0 *  (tau_max/ (w0'*x_active));

NIR_active  = NIR(:, x_ind);


opt.w = w;
opt.num_var = size(NIR_active, 2);
[x_r, S_r, c_r, run_r] = BB_operator(NIR_active, octane, opt);

%reweight2
x_ind1 = find(x_r);
x_active1 = x_r(x_ind1);
w_active1 = w(x_ind1);
w01 = 1./(abs(x_active1) + epsilon);
w1 = w01 *  (tau_max/ (w01'*abs(x_active1)));
%w1 =  w(x_ind1);

NIR_active1  = NIR_active(:, x_ind1);

opt.w = w1;
opt.num_var = size(NIR_active1, 2);
[x_r1, S_r1, c_r1, run_r1] = BB_operator(NIR_active1, octane, opt);

%x_pw = ANEN_helpers(NIR, octane, tau, delta);

x_r_out = zeros(size(x_l1));
x_r_out(x_ind) = x_r;

x_r_out1 = zeros(size(x_l1));
x_r_out1(x_ind(x_ind1)) = x_r1;

figure(25)
clf
hold on
stem(x_l1, '.')
stem(x_r_out, '.')
stem(x_r_out1, '.')
%stem(x_max, '.')
%stem(x_p, '.')
%stem(x_en, '.')
%stem(x_delta, '.')
%stem(x_pw, '.')
%legend('Constrained min (affine)', 'L2 min')
hold off

error_l1 = run_l1.warm_start.bash_manager.get_error(run_l1.warm_start.y);
error_r  = run_r.warm_start.bash_manager.get_error(run_r.warm_start.y);
error_r1  = run_r1.warm_start.bash_manager.get_error(run_r1.warm_start.y);

function Ax = data_A(A, x)
    Ax = A*x;
end

function Atx = data_At(A, x)     
    Atx = A'*x;    
end
