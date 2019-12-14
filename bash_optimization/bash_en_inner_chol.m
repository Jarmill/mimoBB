function [c_bash] = bash_en_inner_chol(c0, tau, K_full, rhs_full, CS)
%performs affine bashing over the EN objective
%Uses low-rank cholesky updates
%
%Input:
%   c0:         Initial set of coefficients on the atoms
%   tau:        Radius of L1-ball/Simplex
%   K_full:     Kernel for full-dimensional affine steps
%   rhs_full:   Right hand side for full affine steps (depends on b)
%   CS:         Whether overall atomic set is centrally symmetric, so that
%               negative values of c are acceptable.
%
%Output:
%   c_bash:         Optimal loading of coefficients over atoms in the set


%create the cholesky manager, and start on interior or exterior as required
CM = chol_manager(K_full, rhs_full, tau, CS);
%y0 = CM.homog(c0);
y0 = c0 / tau;
y = y0;

on_boundary = abs(norm(y0, 1) - 1) < 1e-12; %numerical precision warning?


s = ones(size(y));
s(sign(y) == -1) = -1;


if on_boundary
    CM = CM.full_to_ext(s);
    y = abs(y(1:end-1));
end

%main loop setup
inner_done = 0;
bash_iter = 1;

%signs are handled by CM.update_ext.s

while ~inner_done
    [CM, y, stable_face, lockout, drop_ind] = CM.single_bash(y);
    if stable_face
        if isempty(CM.update_ext)
            %interior, done regardless
            inner_done = 1;
        else
            %exterior, need to check stable-suboptimality
            grad = CM.gradient_full(y);
            s = [CM.update_ext.s; CM.update_ext.s_anchor];            
            if all(sign(grad) == -s)
                %stable-optimal face
                inner_done = 1;
            else
                %stable-suboptimal face
                y =  [CM.update_ext.s .*y; CM.update_ext.s_anchor*(1-sum(y))];
                CM = CM.ext_to_full();
            end
        end
    end
    bash_iter = bash_iter + 1;
end

c_bash = CM.coeff(y);