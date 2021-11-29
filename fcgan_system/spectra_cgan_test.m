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
tau = 100;
%delta = 1.5;

% tau_max = 100;
% 
% opt.num_var = size(NIR, 2);
% opt.tau = tau_max;
% opt.delta = 0;
% opt.norm_type = 1;
% opt.is_complex = 0;
% opt.w = 1;
% opt.export_warm_start = 1;

param = struct;
param.rho=tau;
param.max_nb_atoms=500; % Right now the code assumes that both are the same
param.max_nb_iter=500;
max_iter_fw=500;
param.epsStop=1e-5;
param.debug=false;
param.lmo=@(g) LMO_1d(g, 1, 1);
param.num_var = size(NIR, 2);
param.method = 'asqp';
%BB_1d
% [x_l1, S_l1, c_l1, run_l1] = BB_operator(NIR, octane, opt);
[x, as, hist, iter] = cgan_operator(NIR,octane,param);


figure(25)
clf
hold on
stem(x_l1, '.')
hold off

error_l1 = run_l1.warm_start.bash_manager.get_error(run_l1.warm_start.y);

% function Ax = data_A(A, x)
%     Ax = A*x;
% end
% 
% function Atx = data_At(A, x)     
%     Atx = A'*x;    
% end
