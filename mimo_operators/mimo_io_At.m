function [Atb] = mimo_io_At(b,np, nu, ny, Ns, F, ha, f, W)
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

%legacy (no longer used)
%   Wt: weigthing term in determinant minimization

%b = reshape(b, Ns, ny);


b_time = b(1:(Ns*ny));

if Ns == 0
    Atb = zeros(ny, nu, np);
else
    Atb_time = mimo_At(b_time,np, nu, ny, Ns, F, ha);
    Atb = Atb_time;
    Atb = reshape(Atb, ny, nu, np);
end

b_freq_real = b(Ns*ny + 1:end);
b_freq = complex_fold(b_freq_real, 1);


b_freq = reshape(b_freq, ny, nu, size(f, 1));


%weighting term
%f = Tr(E' W E)
%    W = (E'E)^-{1}

%Atb = zeros(np, nu, ny);
%Atb = reshape(Atb, np, nu, ny);



Atb_freq = mimo_freq_At(b,np, nu, ny, Ns, f, W);

% fwb = zeros(size(Atb));
% 
% %could probably replace with a tensor product
% for j = 1:nu
%     for i = 1:ny
%         %b_curr  = b_freq(:, j, i);
%         %wb = conj(W(:, j, i)).*b_curr;
%         %uwb = conj(U(:, j)) .* wb;
%         %fuwb = f'*uwb;
%         
%         b_curr  = b_freq(i, j, :);
%         wb = squeeze(conj(W(i, j, :)).*b_curr);
%         fwb(i, j, :) = f'* wb;
%         
%         
%         %Atb(i, j, :) = Atb(i, j, :) + fwb;
%         %Atb(:, j, i) = Atb(:, j, i) + fwb;
%         %ind_curr = ny*nu*(0:(np-1))  + nu*(i-1) + j;
%         %Atb(ind_curr) = Atb(ind_curr) + fuwb;
%     end
% end
Atb = Atb + real(Atb_freq);
%Atb = real(Atb);
Atb = squeeze(reshape(Atb, [], 1 ,1));

end