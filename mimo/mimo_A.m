function [Ac] = mimo_A(c,np, nu, ny, Ns, F, ha)
%linear operator for response at output of subsystems with respect to
%input. Used in c -> Ac in error and gradient computation

%make sure this ordering is correct
%c = reshape(c, np, nu, ny);
%reshape doesn't work with sparse matrices


Ac = zeros(Ns, ny);


c0 = [];
rc0 = [];
for i = 1:ny    
    for j = 1:nu
        %c_curr = squeeze(c(:, j, i ));
        ind_offset = (j-1)*np + nu*np*(i-1);
        ind_curr = (1:np) + ind_offset;
        c_curr = c (ind_curr);
        rc_curr = ha*c_curr;
        urc_curr = toeplitzmult2(F{j}, rc_curr);
        Ac(:, i) = Ac(:, i) + urc_curr;
        c0 = [c0; c_curr];
        rc0 = [rc0; rc_curr];
    end    
end

Ac = reshape (Ac, [], 1);



end

