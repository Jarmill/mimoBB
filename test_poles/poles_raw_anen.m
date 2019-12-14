%attempt using an L1 norm fit on the poles
rng(100, 'twister');
Fs = 1;

%p_sys = [0.7];
%p_sys = [0.3+0.5j, 0.3-0.5j, 0.8];
%p_sys = [0.9*exp(2.5j), 0.9*exp(-2.5j)];
%p_sys = [0.8*exp(0.1j), 0.8*exp(-0.1j)];
p_sys = [0.8+(0.1j), 0.8-(0.1j)];
%p_sys = [0.95*exp(0.1j), 0.95*exp(-0.1j)];
%p_sys = [-0.9; 0.5];
%p_sys = [0; 0.9];
%p_sys = [0.65];
%p_sys = [1/sqrt(2)];
%p_sys = [-0.005 + 0.5j; -0.005 - 0.5j];
%p_sys = [0.95j; -0.95j];
%p_sys = [0.7j; -0.7j; 0.2];
%p_sys = [-0.5 + 0.5j, -0.5 - 0.5j];
%p_sys = [-1];
%p_sys = [0.75; 0.95];
%p_sys = [0.4; 0.9];
%p_sys = [0.3; 0.5];
%p_sys = [0.3; 0.7];
%p_sys = [0.8];
%p_sys = [-0.5 + 0.5j, -0.5 - 0.5j, 0.7];

z_sys = [];
k_sys = 1;
sysd = zpk(z_sys, p_sys, k_sys, Fs);

%regularization and output
%tau = 1.8;
%tau = 1.9;      %atomic ball radius
%tau = 2;
%tau = 5;
%tau = 5.5;
%tau = 7;
%tau = 10;
%tau = 15;
tau = 18;
%tau = 20;
%tau = 25;

%delta = 1e-4;   %Elastic Net Regularization
delta = 0;   %Elastic Net Regularization
N = 101;        %time horizon
%y = impulse(sysd, 0:N-1); %output of system
%y = impulse(sysd, 1:N); %output of system

h = impulse(sysd,N-1); 

u = idinput(N,'rbs'); 

%u = zeros(size(h));
%u(1) = 1;

%load u 
%[Wn,zeta] = damp(sys);
%u = sin(Wn(1)*(0:N-1))';
%u = u+ sin(2*Wn(1)*(0:N-1))';
Tu = toeplitz(u,[u(1) zeros(1,N-1)]);
y_true = Tu*h;

noise_level = 0.0;
%noise_level = 0.01;
noise = randn(size(h)) * noise_level;

y = y_true + noise;

%plot(y)
%% define location of poles to be tested
%radius = 5;
%radius = 10;
%radius = 20;
%radius = 30;
radius = 40;
group = 1;



%[A, w, poles] = pole_disk_grid(radius, N, group);
[A_dict, w, poles] = pole_disk_grid(radius, N-1, group);
A = Tu * A_dict;


%% Run BB
use_cvx = 0;
if use_cvx == 1
    cvx_begin    
        M =length(w);
        variable x(M);
        res = A*x - y;
        cost = res'*res;
        norm(w.*x, 1) <= tau;
        
        minimize(cost)
        
    cvx_end
elseif use_cvx == 2
    [x,dual_poly_coeffs,Tu] = ast_cvx(y,tau);
    
else   
    [x, S, c, run_log] = BB_1d(A, y, tau, delta, 1+group, w);    
end

y_rec = A*x;
poles_active = poles(x ~= 0);
c_active = x(x~=0);

%reconstruct the system
spacing_tol = (1/radius)/2;
% sysr = minreal(sysp, spacing/2);

[sysp, sysr] = zpk_from_poles(c_active, poles_active, spacing_tol, Fs);

[g, ng] = gapmetric(sysd, sysp);
[gr, ngr] = gapmetric(sysd, sysr);

figure(2)
clf
hold on
title('System Identification')
plot3(1:N, real(y_rec), imag(y_rec))
plot3(1:N, y, zeros(size(y)), 'k')
hold off
legend('Estimate', 'Ground Truth')