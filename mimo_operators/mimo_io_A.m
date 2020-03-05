function [Ac] = mimo_io_A(c,np, nu, ny, Ns, F, ha, f, W)
%linear operator for response at output of subsystems with respect to
%input. Used in c -> Ac in error and gradient computation

%Frequency domain weighting on individual I/O systems

%Input:
%   c:  coefficient vector to be multiplies
%   [np, nu, ny]:   Number of poles (dict. entries), inputs, and ouputs
%   Ns:             Number of samples
%   F:  FFT for Toe plitz of input u
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

Ac_time = mimo_A(c, np, nu, ny, Ns, F, ha);


Nf = size(W, 3);
%Ac_freq = zeros(Nf, nu, ny);
Ac_freq = zeros(ny, nu, Nf);

%c = reshape(full(c), [np, nu, ny]);
c = reshape(full(c), [ny, nu, np]);

%no addition going on here
for i = 1:ny
    for j = 1:nu
        %I think this is where I was indexing incorrectly 
%        ind_curr = ny*nu*(0:(np-1)) + nu*(i-1) + j;

        
%        c_curr = c (ind_curr);
        %c_curr = c(:, j, i);
        c_curr = squeeze(c(i, j, :));
        w_curr = squeeze(W(i, j, :));
        fc = f * c_curr;
        
        %ufc = U(:, j) .* fc;
        %wfc = W(:, j, i).*ufc;
        %wfc = W(:, j, i) .* fc;
        %Ac_freq(:, j, i) = wfc;
        wfc = w_curr .* fc;
        Ac_freq(i, j, :) = wfc;
    end
end

Ac_freq =  squeeze(reshape(Ac_freq, [], 1, 1));
Ac_freq_real = complex_unfold(Ac_freq, 1);

Ac = [Ac_time; Ac_freq_real];