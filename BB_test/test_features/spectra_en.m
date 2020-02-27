load 'spectra'

tau = 100;
delta = 1;

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
Af  = @(x) data_A(NIR, x);
Atf = @(x) data_At(NIR, x);

[x_l1, S_l1, c_l1] = BB_operator(Af, Atf, octane, opt);
%[x_max, S_max, c_max] = BB_1d(NIR, octane, tau_max, delta, Inf);
%[x_p, S_p, c_p] = BB_1d(NIR, octane, tau_p, delta, p);


%x_pw = ANEN_helpers(NIR, octane, tau, delta);

figure(25)
clf
hold on
stem(x_l1, '.')
%stem(x_max, '.')
%stem(x_p, '.')
%stem(x_en, '.')
%stem(x_delta, '.')
%stem(x_pw, '.')
%legend('Constrained min (affine)', 'L2 min')
hold off


function Ax = data_A(A, x)
    Ax = A*x;
end

function Atx = data_At(A, x)     
    %if nargin == 3          
    %    Atx =A(:, ind)'*x;
    %else
        Atx = A'*x;
    %end        
end
