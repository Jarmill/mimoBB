function [Atb] = mimo_output_At(b,np, nu, ny, Ns, F, ha, f, U, W, Wt)
%adjoint linear operator for response at output of subsystems with respect 
%to input. Used in b -> A'b in gradient computation

%Input:
%   c:  coefficient vector to be multiplies
%   [np, nu, ny]:   Number of poles (dict. entries), inputs, and ouputs
%   Ns:             Number of samples
%   F:  FFT for Toeplitz of input u
%   ha: Pole dictionary in time
%   f:  Pole dictionary in frequency
%   W:  Weighting functions for each I/O pair
%   U:  FFT of input u
%   Wt: weigthing term in determinant minimization

%b = reshape(b, Ns, ny);


b_time = b(1:(Ns*ny));


if nargin >= 11
    b_time = b_time * Wt;
end

b_freq_real = b(Ns*ny + 1:end);
b_freq = complex_fold(b_freq_real, 1);


b_freq = reshape(b_freq, Ns, ny);

Atb_time = mimo_At(b_time,np, nu, ny, Ns, F, ha);
Atb = Atb_time;

%weighting term
%f = Tr(E' W E)
%    W = (E'E)^-{1}
for i = 1:ny
    for j = 1:nu    
        b_curr  = b_freq(:, i);
        wb = conj(W(:, j)).*b_curr;
        uwb = conj(U(:, j)) .* wb;
        fuwb = f'*uwb;
        
        %Atb(:, j, i) = rub_curr;
        ind_curr = ny*nu*(0:(np-1))  + nu*(i-1) + j;
        Atb(ind_curr) = Atb(ind_curr) + fuwb;    
    end
end

Atb = real(Atb);
%Atb = squeeze(reshape(Atb, [], 1 ,1));

end