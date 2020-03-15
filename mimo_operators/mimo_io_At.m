function [Atb] = mimo_io_At(b,np, nu, ny, Ns, Tu, ha, fa, W, Wt)
%adjoint linear operator for response at output of subsystems with respect 
%to input. Used in b -> A'b in gradient computation

%Input:
%   c:  coefficient vector to be multiplies
%   [np, nu, ny]:   Number of poles (dict. entries), inputs, and ouputs
%   Ns:             Number of samples
%   Tu:  Toeplitz of input u
%   ha: Pole dictionary in time
%   fa:  Pole dictionary in frequency
%   W:  Weighting functions for each I/O pair
%   Wt: weigthing term in determinant minimization

%b = reshape(b, Ns, ny);


b_time = b(1:(Ns*ny));


if nargin >= 10
    b_time = b_time * Wt;
end

Atb_time = mimo_At(b_time, np, nu, ny, Ns, Tu, ha);
Atb = Atb_time;

b_freq_real = b(Ns*ny + 1:end);
b_freq = complex_fold(b_freq_real, 1);
b_freq = reshape(b_freq, ny, nu, size(fa, 1));

%weighting term
%f = Tr(E' W E)
%    W = (E'E)^-{1}

%Atb = zeros(np, nu, ny);
%Atb = reshape(Atb, np, nu, ny);
Atb = reshape(Atb, ny, nu, np);
fwb = zeros(size(Atb));

%could probably replace with a tensor product
for j = 1:nu
    for i = 1:ny
        %b_curr  = b_freq(:, j, i);
        %wb = conj(W(:, j, i)).*b_curr;
        %uwb = conj(U(:, j)) .* wb;
        %fuwb = f'*uwb;
        
        b_curr  = b_freq(i, j, :);
        wb = squeeze(conj(W(i, j, :)).*b_curr);
        fwb(i, j, :) = fa'* wb;
        
        
        %Atb(i, j, :) = Atb(i, j, :) + fwb;
        %Atb(:, j, i) = Atb(:, j, i) + fwb;
        %ind_curr = ny*nu*(0:(np-1))  + nu*(i-1) + j;
        %Atb(ind_curr) = Atb(ind_curr) + fuwb;
    end
end
Atb = Atb + real(fwb);
%Atb = real(Atb);
Atb = squeeze(reshape(Atb, [], 1 ,1));

end