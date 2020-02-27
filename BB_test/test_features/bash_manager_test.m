rng(30)
N = 7;

A = rand(10, N);
b = 10*rand(10, 1);
tau = 10;
%tau = 1.5;
delta = 1e-3;
norm_type = 1;

BM = bash_manager(b, tau, delta, norm_type);

I = eye(N);
S1  = I(:, 1:2);
AS1 = A*S1;

S2  = I(:, 3:5);
AS2 = A*S2;

S3 = I(:, 6:7);
AS3 = A*S3;

BM = BM.add_atoms(S1, AS1);
BM = BM.add_atoms(S2, AS2);

s = [1; 1; -1; 1; -1];
BM = BM.sign_switch(s);

% BM = BM.full_to_ext();
%
%y = BM.solve_system();
BM = BM.delete_index(2);
%y = BM.solve_system();

BM = BM.full_to_ext();
%BM = BM.ext_to_full();


BM = BM.add_atoms(S3, AS3);
% 
% BM = BM.delete_anchor();

BM = BM.ext_to_full();


% BM = BM.add_atoms(S1, AS1);
% BM = BM.sign_switch([1; -1]);
% BM = BM.delete_index(1);
% BM = BM.full_to_ext();
% BM = BM.add_atoms(S2, AS2);
% BM = BM.delete_anchor();
% BM = BM.add_atoms(S3, AS3);

% [BM, y_bash] = BM.bash([], I, A*I);
% x_bash = BM.get_x(y_bash);
% grad_bash = BM.grad(A, y_bash)