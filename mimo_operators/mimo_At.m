function [Atb] = mimo_At(b,np, nu, ny, Ns, F, ha, Wt)
%adjoint linear operator for response at output of subsystems with respect 
%to input. Used in b -> A'b in gradient computation

b = reshape(b, Ns, ny);
Atb = zeros( ny , nu, np );
%Atb = zeros(np*nu*ny, 1);

%weighting term
%f = Tr(E' W E)
%    W = (E'E)^-{1}
if nargin >= 8
    b = b * Wt;
end

for j = 1:nu
    for i = 1:ny
        b_curr  = b(:, i);
        ub_curr = toeplitzmult2(conj(F{j}), b_curr);
        rub_curr = ha'*ub_curr;
        %Atb(:, j, i) = rub_curr;
        Atb(i, j, :) = rub_curr;
        %ind_curr = ny*nu*(0:(np-1))  + nu*(i-1) + j;
        %Atb(ind_curr) = rub_curr;
    end
end

Atb = squeeze(reshape(Atb, [], 1 ,1));

end