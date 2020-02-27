function [x_out, valid, alpha_max] = intersect_l1_ball(s, x1, x2, tau)
%INTERSECT_L1_BALL
%Given a point x1 on the L1 ball and x2 possibly outside it, this routine
%finds the point on the path x1 -> x2 as close to x2 as possible but still
%on/within the L1-ball.
%
%Numerical error causes issues when x1 and x2 are on the same extended
%face. Split into two parts

%s = sign(x1);
on_boundary = (abs(norm(x1, 1)-tau) < 1e-8);
same_face = (abs(s'*x2 - tau) < 1e-8);

if on_boundary && same_face
    %Exterior step being conducted, x1 and x2 are on the same extended
    %face. Boundary processing can use the simplex intersection routine.
    x1_pos = s.*x1;
    x2_pos = s.*x2;
    [x_out_pos, valid, alpha_max] = intersect_simplex(x1_pos, x2_pos, tau);
    %x_out = x1 + alpha_max*(x2 - x1);
    x_out = s.*x_out_pos;
else
    [x_out, valid, alpha_max] = intersect_l1(x1, x2, tau);
end