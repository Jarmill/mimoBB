function [Atb] = mimo_sys_At(b,np, nu, ny, Ns, F, ha, f, W, Wt)
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
%   Wt: weigthing term in determinant minimization

b = reshape(b, Ns, ny);

b_time = b(1:(Ns*ny));

b_freq = (Ns*ny + 1:end);

b_freq = reshape(b_freq, Ns, nu, ny);

Atb_time = mimo_sys_At(b_time,np, nu, ny, Ns, F, ha, Wt);
Atb = Atb_time;

%weighting term
%f = Tr(E' W E)
%    W = (E'E)^-{1}

for j = 1:nu
    for i = 1:ny
        b_curr  = b_freq(:, j, i);
        wb = W(:, j, i)'.*b_curr;
        fwb = f'*wb;
        
        %Atb(:, j, i) = rub_curr;
        ind_curr = ny*nu*(0:(np-1))  + nu*(i-1) + j;
        Atb(ind_curr) = Atb(ind_curr) + fwb;
    end
end

%Atb = squeeze(reshape(Atb, [], 1 ,1));

end