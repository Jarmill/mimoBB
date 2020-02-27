p = [0.3; -0.8; 0.3 + 0.6j; 0.7j; 0.9 + 0.1j; 0.6 + 0.8j];
N = 20;

[A, w] = pole_matrix_upper(p, N);