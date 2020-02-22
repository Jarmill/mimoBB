function [M_out] = complex_fold(M_in, dims)
%COMPLEX_FOLD fold real matrix into matrix of complex numbers
%   dims = 1: complex numbers in columns [real; imag]
%   dims = 2: complex numbers in 2x2 matrices [real -imag; imag real]

if nargin < 2
    dims = 1;
end

if dims == 1
    M_out = M_in(1:2:end, :) + 1j * M_in(2:2:end, :);
else
    M_out = M_in(1:2:end, 1:2:end) + 1j * M_in(2:2:end, 1:2:end);
end

end