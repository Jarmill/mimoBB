function [Atb_freq] = mimo_freq_At(b,np, nu, ny, Ns, f, W)
%adjoint linear operator for response at output of subsystems with respect 
%to input. Used in b -> A'b in gradient computation
%only the frequency term

%Input:
%   c:  coefficient vector to be multiplies
%   [np, nu, ny]:   Number of poles (dict. entries), inputs, and ouputs
%   Ns:             Number of samples
%   W:  Weighting functions for each I/O pair

%figure out this indexing later
b_freq_real = b(Ns*ny + 1:end);
b_freq = complex_fold(b_freq_real, 1);

b_freq = reshape(b_freq, ny, nu, size(f, 1));


%weighting term
%f = Tr(E' W E)
%    W = (E'E)^-{1}

Atb_freq = zeros(np, nu, ny);
%Atb = reshape(Atb, np, nu, ny);
Atb_freq = reshape(Atb_freq, ny, nu, np);
% fwb = zeros(size(Atb));

%could probably replace with a tensor product
for j = 1:nu
    for i = 1:ny
        
        b_curr  = b_freq(i, j, :);
        wb = squeeze(conj(W(i, j, :)).*b_curr);
        Atb_freq(i, j, :) = f'* wb;
        

    end
end

% Atb = squeeze(reshape(Atb, [], 1 ,1));

end