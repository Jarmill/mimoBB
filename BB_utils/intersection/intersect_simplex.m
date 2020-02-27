function [ x_out, valid, alpha_max ] = intersect_simplex(x1, x2, tau)
%INTERSECT_simplex finds the a point inside/on the simplex on the path between
%x1 and x2, as close to x2 as possible.
%
%If the path between x1 and x2 goes through the simplex, then valid = 1
%Otherwise, valid = 0
%
%Formerly known as "the stupid algorithm". Still the stupid algorithm.
%
%x2 has a lower error (function value) than x1, so by strong convexity, all
%points on path between x1->x2 has a lower error than f(x1) (as long as
%f(x2) < f(x1), Jensen's inequality in action). 
%
%direction of motion and stepsize (alpha)
%linear interpolation x_out = x1 + alpha*(x2-x1)
%
% Assumes for now that x1 is a feasible point (within or on the simplex)
w = x2 - x1;
x_out = x2;
%s = sign(x_out);
s = ones(size(x_out));

%intersection with sidewall
if s'*x_out - tau <= 1e-9 && all(x_out >= 0)
    alpha_max = 1;
    valid = 1;
else
    %sidewall intersection
    if any(w(x1 == 0) < 0)
        x_out_sidewall = x1;
        alpha_max_sidewall = 0;
        %valid = 0;
    else
        alpha_list = -x1./w;
        alpha_max_sidewall = min(alpha_list(alpha_list > 0));
        if isempty(alpha_max_sidewall)
           alpha_max_sidewall = 1;
        end
        %valid = 1;
        %x_out_sidewall = x1 + alpha_max_sidewall*w;
    end

    %intersection with face
    on_boundary = abs(tau - s'*x1) < 1e-9;
    if on_boundary
        alpha_max = alpha_max_sidewall;
    else
        alpha_max_face = (tau - s'*x1)/(s'*w);    
        
        if (alpha_max_face < alpha_max_sidewall) && (alpha_max_face >= 0)
            alpha_max = alpha_max_face;
        else
            alpha_max = alpha_max_sidewall;
        end
        %this min call doesn't
        %alpha_max = min(alpha_max_sidewall, alpha_max_face);
    end
    valid = (alpha_max > 0);
    %x_out_face = x1 + alpha_max_face*w;    
    %norm_gap = norm(x_out, 1) - tau;
    
    x_travel = x1 + alpha_max*w;
    
    %lockout on the anchor, no other way to detect it (I think)
    if on_boundary && (s'*w > 0)
        x_out = x1;
        alpha_max = 0;
    else
        x_out = x_travel;
    end
end


end