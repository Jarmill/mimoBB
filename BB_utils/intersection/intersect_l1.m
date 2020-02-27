function [ x_out, valid, alpha_max ] = intersect_l1(x1, x2, tau)
%INTERSECT_L1 finds the a point inside/on the L1-ball on the path between
%x1 and x2, as close to x2 as possible.
%
%If the path between x1 and x2 goes through the L1 ball, then valid = 1
%Otherwise, valid = 0
%
%Formerly known as "the stupid algorithm"
%
%x2 has a lower error (function value) than x1, so by strong convexity, all
%points on path between x1->x2 has a lower error than f(x1) (as long as
%f(x2) < f(x1), Jensen's inequality in action). 

%floating point/numerical error sucks.
%Write a new intersection routine for L1.
%x1 is feasaible, x2 may not be. If x2 is on same extended face as x1 (and
%x1 on boundary) do a simplex/like line search to find the point of
%intersection. In all other circumstances (x2 inside, x1 inside, x1 on
%boundary but x2 does not share extended face), use the stupid algorithm to
%find the point of intersection


%direction of motion and stepsize (alpha)
%linear interpolation x_out = x1 + alpha*(x2-x1)
w = x2 - x1;
alpha_max = 1;

%output point x_out. Try out x2 first. 
%x_out = x1 + alpha_max*w;
x_out = x2;
norm_gap= norm(x_out, 1) - tau;

while 1
    %check if the output point is valid (on or inside L1 ball)
    if norm_gap < 1e-8 
        valid = 1;
        break;
    end
    
    alpha_max_old = alpha_max;
    
    %Otherwise, find a point intersecting a plane splayed from the L1-ball    
    s = sign(x_out);
    alpha_max = (tau - s'*x1)/(s'*w);
    x_out = x1 + alpha_max*w;    
    norm_gap = norm(x_out, 1) - tau;
    
    %plot(x_out(1), x_out(2), 'xk')
    
    %This is a sign that the path does not intersect the L1-ball
    %Either the alpha goes up (oscillating between 2+ infeasible points,
    %with on different orthants, or x1->x2 on same orthant points towards 
    %L1-ball but does not reach it), or alpha < 0 (x1->x2 on same orthant
    %but x2 is further away than x1).
    if alpha_max > alpha_max_old || alpha_max < 0
        valid = 0;
        break;
    end        
end



end