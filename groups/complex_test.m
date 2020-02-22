rng(40);

N = 5;
delta = 0.01;
%tau = 10;
%tau = 3;
tau = 2;

Ar = randn(N);
Ai = randn(N);

br = 5*randn(N, 1);
bi = 5*randn(N, 1);


A = Ar + 1j*Ai;
b = br + 1j*bi;

% K = A'*A + delta*eye(N);
% 
% x = A \ b;
% x_tikh = K \ (A'*b);


%atom testing ground
% xn = x_tikh./abs(x_tikh);
% S = diag(xn);
% AS = A*S;
% 
% KS = (AS'*AS) + delta*(S'*S);
% rhs = AS'*b/tau;
% 
% c = KS \ rhs;

%S = bag_l1(-A'*b, zeros(size(x)), 3, 1);
% S = bag_complex(-A'*b, zeros(size(x)), 3, tau, 1);
% AS = A*S;
% KS = (AS'*AS) + delta*(S'*S);
% rhs = AS'*b/tau;
% 
% c = KS \ rhs;
% 
% AS_w = complex_unfold(AS, 1);
% S_w = complex_unfold(S, 1);
% S_w2 = complex_unfold(S, 2);
% b_w = complex_unfold(b, 1);
% 
% K_w =  (AS_w'*AS_w) + delta*(S_w'*S_w);
% rhs_w = AS_w'*b_w/tau;
% 
% c_w = K_w \ rhs_w
% 
% x_w = S*c_w

% theta = 2*pi*rand(N, 1);
% St = diag(exp(theta*1j));
% ASt = A*St;
% KSt = (ASt'*ASt) + delta*(St'*St);
% rhst = ASt'*b/tau;
% 
% ct = KSt \ rhst;
% ctr = real(ct);
% 
% xc = tau * S*c;
% xt = tau * St*ct;
% xtr = tau* St*ctr;



%x = zeros(size(z5));

%[S_bag, N_added] = bag_l1(z5, x, 3);
%[x_bb, S_bb, c_bb] = BB_1d(A, b, tau, delta, 1);
%[x_bb, S_bb, c_bb] = BB_1d(A, b, tau, delta, Inf);
[x_bb, S_bb, c_bb] = BB_1d(A, b, tau, delta, 1.5);


grad_bb = A'*(A*x_bb - b) + delta*x_bb;

phase_diff = angle(x_bb) - angle(-grad_bb);