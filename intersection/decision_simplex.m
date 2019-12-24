function [x_out, outside_flag, drop_ind] = decision_simplex(x, x_aff, tau, on_boundary)

%makes a decision on the new iterate x_out given x and x_aff (f(x) >
%f(x_aff)) over the simplex.
%only supports naive affine steps for now

%outside_flag
%   0: on interior, feasible
%   1: on boundary of L1-ball, feasible
%   2: outside boundary on affine subspace, infeasible

if on_boundary
    %exterior things here
    %copied from exterior_step
    [x_out, ~, alpha_max] = intersect_simplex(x, x_aff, tau);

    if alpha_max > 1e-9 && alpha_max < 1
        %unstable face
        %more inefficient than previous, but numerically accurate?
        drop_ind = find(x_out ~= 0 & abs(x_out) < 1e-8, 1);
        [~, i] = min(abs(x_out(drop_ind)));
        x_out(drop_ind(i)) = 0;
        outside_flag = 2;
    else

        if alpha_max < 1e-9
            outside_flag = 0; %invalid affine step, stays in same position. use fewer dimensions
        elseif alpha_max == 1
            outside_flag = 1; %valid exterior affine step
        else
            outside_flag = 2; %outside region
        end
        drop_ind = -1;
    end
    %outside_flag = (alpha_max ~= 1) + 1;
else
    %interior things
    [x_out, ~, alpha_max] = intersect_simplex(x, x_aff, tau);
    outside_flag = (alpha_max > 0)*(1 + (alpha_max < 1)); %awful conditional logic
    drop_ind = -1;
end