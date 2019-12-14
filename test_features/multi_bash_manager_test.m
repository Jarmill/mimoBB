rng(30)
N = 7;

A = rand(10, N);
b = 10*rand(10, 1);
tau = 10;
%tau = 1.5;
delta = 1e-3;
norm_type = 1;


%constraints
con1 = struct;
con1.norm_type = 'simp';
con1.tau = 2;
con1.w = 1;
con1.index = [1, 2, 3, 4];

con2 = struct;
con2.norm_type = 'simp';
con2.tau = 4;
con2.w = 1;
con2.index = [5, 6, 7];


cons = {con1, con2};
FCFW = 1;

% function testing
[x_final, S_final, c_final, atom_source_final, run_log] = BB_1d_multi(A, b, cons, delta);


% %manual testing
% BM = multi_bash_manager(b, cons, delta, FCFW);
% 
% %start testing
% x = zeros(N, 1);
% Ax = zeros(10, 1);
% c0 = [];
% grad = A'*(Ax - b) + delta*x;
% 
% %iter 1
% [BM, S_bag, atom_source_bag, DG] = BM.bag_atoms(grad, x, 1);
% 
% AS_bag = A*S_bag;
% 
% [BM, c1] = BM.bash(c0, S_bag, AS_bag, atom_source_bag);
% 
% %iter 2
% x = BM.get_x(c1);
% Ax = BM.get_Ax(c1);
% grad = BM.grad(A, c1);
% [BM, S_bag, atom_source_bag, DG] = BM.bag_atoms(grad, x, 1);
% 
% AS_bag = A*S_bag;
% [BM, c2] = BM.bash(c1, S_bag, AS_bag, atom_source_bag);
% 
