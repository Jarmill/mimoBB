function [ A, w, poles] = pole_matrix_upper( p, N, group)
%Vandermonde matrix on the poles
%performs task of finding impulse responses components to horizon
%only use poles on the upper half plane

%p: set of poles in complex plane (must include conjugates)
%N: horizon, number of samples
%group:     0: Standard L1 fit
%           1: Group Lasso

%separate out exponential and trigonometric poles
p_exp = p(imag(p) == 0);
p_trig = p(imag(p) ~= 0);

k_exp = length(p_exp);
k_trig = length(p_trig);
k_all = k_exp + 2*k_trig;
%in total there are k_exp + 2*k_trig atoms

%response of 1/(z-p) = r^n (cos(n theta) + j sin(n theta))
r_trig = abs(p_trig);
theta = angle(p_trig)';

%raise r to powers
%index_N = ones(1,N-1-zero_pad);
h_exp = (p_exp.^(0:N-1))';

mag = (r_trig.^(0:N-1))';
angles = (0:N-1)'*theta; %n*theta

%optional step
angles = mod(angles, 2*pi);


%trig_angles = zeros(size(angles));
%trig_angles(:, imag(p) >= 0) = cos(angles(:, imag(p) >= 0));
%trig_angles(:, imag(p) < 0) = sin(angles(:, imag(p) < 0));

trig_angles_cos = cos(angles);
trig_angles_sin = sin(angles);

mag(abs(mag) < 1e-15) = 0;
trig_angles_cos(abs(trig_angles_cos) < 1e-15) = 0;
trig_angles_sin(abs(trig_angles_sin) < 1e-15) = 0;

%assemble into vandermonde matrix
vandermonde_cos = mag .* trig_angles_cos;
vandermonde_sin = mag .* trig_angles_sin;
vandermonde_exp = h_exp;

A = [zeros(1, k_all); vandermonde_cos vandermonde_sin vandermonde_exp];
%A = [vandermonde_cos vandermonde_sin vandermonde_exp];

%now find the group lasso weights
delta = 1e-6;

%w_exp = ones(size(p_exp));
%w_trig = ones(size(p_trig));

w_exp_0  = 1./(1 - abs(p_exp).^2 + delta);
w_trig_0 = 1./(1 - abs(p_trig).^2 + delta);

w_exp = (1 - abs(p_exp).^(2*N - 2))./(1 - abs(p_exp).^2);
w_trig = (1 - abs(p_trig).^(2*N - 2))./(1 - abs(p_trig).^2);

w_exp(1-abs(p_exp) <= 1e-8) = 1/delta;
w_trig(1-abs(p_trig) <= 1e-8) = 1/delta;

poles = [p_trig; conj(p_trig); p_exp];
if group
    w = struct;
    w.groups = cell(k_trig + k_exp, 1);
    w.pole_groups = cell(k_trig + k_exp, 1);
    w.weights = [w_trig; w_exp];
    w.k_trig = k_trig;
    w.k_exp = k_exp;    

    for i = 1:k_trig
        w.groups{i} = [i, i + k_trig]; 
        w.pole_groups{i} = [poles(i), poles(i + k_trig)]; 
    end

    %index magic
    for i = (k_trig + 1):(k_trig + k_exp)
        w.groups{i} = k_trig + i;
        w.pole_groups{i} = poles(k_trig + i);
    end
else
    w = [w_trig; w_trig; w_exp];
end

end

