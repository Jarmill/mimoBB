function [Ac] = mimo_sys_A(c,np, nu, ny, Ns, F, ha, f, W)
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
%   W:  Weighting functions for each I/O pair


%Ac: [time domain
%     freq response 11
%     freq response 12 
%     ...]

%make sure this ordering is correct
%c = reshape(c, np, nu, ny);
%reshape doesn't work with sparse matrices

%different order of variables, should be more amenable to randomization

Ac_time = mimo_A2(c, np, nu, ny, Ns, F, ha);

Ac_freq = zeros(Ns, nu, ny);

%no addition going on here
for i = 1:ny
    for j = 1:nu
        ind_curr = ny*nu*(0:(np-1)) + nu*(i-1) + j;
        
        c_curr = c (ind_curr);
        fc = f * c_curr;
        wfc = W(:, j, i).*fc;
        Ac_freq(:, j, i) = wfc;
    end
end

Ac = [Ac_time; squeeze(reshape(Ac_freq, [], 1, 1))];