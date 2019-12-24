function [x_out, step_code, drop_ind] = decision_l1(s, x, x_aff, tau, on_boundary)

%makes a decision on the new iterate x_out given x and x_aff (f(x) >
%f(x_aff)) over the L1-ball
%only supports naive affine steps for now

%step_code
%on interior
%   0: interrupted by boundary
%   1: stable on the interior
%on exterior
%   0: lockout (infeasible, cannot fix)
%   1: stable face (feasible)
%   2: interrupted by a boundary (infeasible, keep bashing)



    if on_boundary
        %exterior things here
        %copied from exterior_step
        %[x_out, ~, alpha_max] = intersect_l1(x, x_aff, tau);
        [x_out, ~, alpha_max] = intersect_l1_ball(s, x, x_aff, tau);
        
        if alpha_max > 0  && alpha_max < 1
            %if alpha_max > 1e-9  && alpha_max < 1
            %unstable face
            %more inefficient than previous, but numerically accurate?
            drop_ind = find(x_out ~= 0 & abs(x_out) < 1e-8, 1);
            [~, i] = min(abs(x_out(drop_ind)));
            x_out(drop_ind(i)) = 0;
            step_code = 2;
        else

            if alpha_max <=  1e-9 %numerical error may result in alpha = -1.4e-16
                step_code = 0; %LOCKOUT: x stays in same position.
            elseif alpha_max == 1
                step_code = 1; %valid exterior affine step
            else
                step_code = 2; %outside region
            end
            drop_ind = -1;
        end
        %step_code = (alpha_max ~= 1) + 1;
    else
        %interior things
        [x_out, ~, alpha_max] = intersect_l1(x, x_aff, tau);
        step_code = (alpha_max < 1);
        drop_ind = -1;
    end
% end