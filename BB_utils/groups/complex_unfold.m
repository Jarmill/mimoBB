function [M_out] = complex_unfold(M_in, dims)
%COMPLEX_UNFOLD Unfold matrix of complex numbers into real matrices
%   dims = 1: complex numbers in columns [real; imag]
%   dims = 2: complex numbers in 2x2 matrices [real -imag; imag real]

if nargin < 2
    dims = 1;
end

if dims == 1
    if issparse(M_in)
        M_out = sparse(size(M_in, 1)*2, size(M_in, 2));
    else
        M_out = zeros(size(M_in, 1)*2, size(M_in, 2));
    end
    M_out(1:2:end, :) = real(M_in);
    M_out(2:2:end, :) = imag(M_in);
else
    if issparse(M_in)
        M_out = sparse(size(M_in, 1)*2, size(M_in, 2)*2);
    else
        M_out = zeros(size(M_in, 1)*2, size(M_in, 2)*2);
    end
    M_out(1:2:end, 1:2:end) = real(M_in);
    M_out(2:2:end, 1:2:end) = imag(M_in);
    M_out(1:2:end, 2:2:end) = -imag(M_in);
    M_out(2:2:end, 2:2:end) = real(M_in);
end



end

