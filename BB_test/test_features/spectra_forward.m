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

%BB_1d
[x_l1, run_l1] = BB_forward(NIR, octane, opt);
