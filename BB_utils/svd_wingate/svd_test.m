N = 10;
d = 3;

%try adding a column to A
A = randi(10, N, d);
[U, S, V] = svd(A, 'econ');


c = randi(10, N, 1);

%add column
% i = d+1;
% b = zeros(d+1, 1);
% b(i) = 1;
% a = c;
% Vz = [V; zeros(1, d)];
% [Up, Sp, Vp] = rank_one_svd_update(U, S, Vz, a, b, 0);
% 
% 
% A2 = Up*Sp*Vp';
% A2_true = [A c];
% 
% norm(A2 - A2_true)

%delete column
i = 2;
a = -A(:, i);
b = zeros(d, 1);
b(i) = 1;
[Up, Sp, Vp] = rank_one_svd_update(U, S, V, a, b, 0);

Upm = Up(:, 1:end-1);
Spm = Sp(1:end-1, 1:end-1);
Vpm = Vp(:, 1:end-1);
