function [ en_err] = en_abs( A, x, b, delta )
%EN_ABS absolute elastic net error
%   used to find elastic net error of a point x

%err = norm(A*x - b)^2 /2 + delta/2 * norm(x)^2;

%least squares and l2 penalty
lsq = sum((A*x - b).^2, 1);
l2  = sum(x.^2, 1);

%weighted sum of the two errors
en_err = 0.5*(lsq + delta * l2);

end

