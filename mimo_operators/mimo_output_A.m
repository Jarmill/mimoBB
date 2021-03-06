function [Ac] = mimo_output_A(c,np, nu, ny, Ns, F, ha, f, U, W)
%linear operator for response at output of subsystems with respect to
%input. Used in c -> Ac in error and gradient computation

%Frequency domain weighting on individual I/O systems

%Input:
%   c:  coefficient vector to be multiplies
%   [np, nu, ny]:   Number of poles (dict. entries), inputs, and ouputs
%   Ns:             Number of samples
%   F:  FFT for Toeplitz of input u
%   ha: Pole dictionary in time
%   f:  Pole dictionary in frequency
%   U:  FFT of input u
%   W:  Weighting functions for each I/O pair


%Ac: [time domain
%     freq response 11
%     freq response 12 
%     ...]

%make sure this ordering is correct
%c = reshape(c, np, nu, ny);
%reshape doesn't work with sparse matrices

%different order of variables, should be more amenable to randomization

Ac_time = mimo_A(c, np, nu, ny, Ns, F, ha);

Ac_freq = zeros(Ns, ny);

%no addition going on here
for i = 1:ny
    Ac_freq_curr = zeros(Ns, 1);
    for j = 1:nu
        ind_curr = ny*nu*(0:(np-1)) + nu*(i-1) + j;
        
        c_curr = c (ind_curr);
        fc = f * c_curr;
        ufc = U(:, j) .* fc;
        wfc = W(:, i).*ufc;
        Ac_freq_curr = Ac_freq(:, i) + wfc;
    end
    Ac_freq(:, i) = Ac_freq_curr;
end

Ac_freq =  reshape(Ac_freq, [], 1);
Ac_freq_real = complex_unfold(Ac_freq, 1);

Ac = [Ac_time; Ac_freq_real];