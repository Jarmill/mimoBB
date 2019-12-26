function [Atb] = mimo_At(b,np, nu, ny, Ns, F, ha)
%adjoint linear operator for response at output of subsystems with respect 
%to input. Used in b -> A'b in gradient computation

b = reshape(b, Ns, ny);
Atb = zeros( np, nu, ny);


for j = 1:nu
    for i = 1:ny
        b_curr  = b(:, i);
        ub_curr = toeplitzmult2(conj(F{j}), b_curr);
        rub_curr = ha'*ub_curr;
        Atb(:, j, i) = rub_curr;
    end
end

Atb = squeeze(reshape(Atb, [], 1 ,1));

end