%min 1/2 ||Ax-b||_2^2 such that ||x||_p <= tau
%where p is either 1 or infinity


A = eye(2);
b = [3; 0.5];
tau = 1;

H = A'*A;
f = -A'*b;

%OLS

x_ols = A \ b;

%OLS from sover
opts = mpcActiveSetOptions;
[x_mpc_ols] = mpcActiveSetSolver(H, f, zeros(0, 2), zeros(0, 1), zeros(0, 2), zeros(0, 1), false(0, 1), opts)

Ainf_hard = kron(eye(2), [-1;1]);
binf_hard = ones(4,1);
[x_mpc_linf_hard] = mpcActiveSetSolver(H, f, Ainf_hard, binf_hard, zeros(0, 2), zeros(0, 1), false(4, 1), opts)

Hlift = blkdiag(H, 0);
flift = [f; 0];
Ainf_lift = [kron(eye(2), [-1;1]), ones(4, 1); 0,0,1];
binf_lift = [zeros(4,1); tau];
[x_mpc_linf_lift] = mpcActiveSetSolver(Hlift, flift, Ainf_lift, binf_lift, zeros(0, 3), zeros(0, 3), false(5, 1), opts)
