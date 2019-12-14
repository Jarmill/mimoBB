function [c_bash, lockout, lockout_dims, on_boundary_out] = bash_en_inner(c0, tau, K_full, rhs_full, CS)
%performs affine bashing over the EN objective
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
%   lockout:        If the set has locked out, meaning that too many atoms 
%                   have been added from the bag.
%   lockout_dims:   Newly added dimensions that got locked out

on_boundary = abs(norm(c0, 1) - tau) < 1e-9; %numerical precision warning?
inner_done = 0;
lockout = 0;
c = c0;

%signs of the face of the L1-ball
%all atoms start out positive coming from the bag.
s = sign(c);
s(s == 0) = 1;

first_iter = 1; %first iteration

take_interior = ~on_boundary;

%cholesky decomposition. Weird affine magic.
USE_CHOL = 0;
if USE_CHOL
    L_int = chol(K_full);
end

USE_SLIDE = 0;

while ~lockout && ~inner_done
    %loop through affine steps until a stable-optmial face is hit
    
    if take_interior
        %interior affine steps
        %either starting on interior, or on exterior on a stable-suboptimal
        %face. Can never start on a stable-suboptimal face though.
        
        if first_iter
            %on interior, and making decision
            c_full = K_full \ rhs_full;
            c_full(abs(c_full) < 1e-12) = 0; %numerical cleanup (if needed, can disable)
        else
            %handles stable-suboptimal faces
            [B, ~] = affine_ext_representation(s, 1, 0);
            
            if USE_CHOL
                %rhs_full_curr = B'*rhs_full;
                %y_full = affine_chol_solve(B, L_int, rhs_full_curr);
                c_full = affine_chol_solve(B, p, K_full, L_int, rhs_full);
            else
                
                %if K_full is not full rank, this may cause trouble.
                rhs_full_curr = B'*rhs_full / tau;
                K_full_curr = B'*K_full*B;
                y_full = K_full_curr \ rhs_full_curr;
                c_full = tau* B*y_full;
            end                             
            %y_full = affine_chol_solve(B, L_int, rhs_full); 
            
        end
        
        
        if CS
            [c_out, valid_affine_step] = decision_l1(s, c, c_full, tau, 0);
            %lockout is impossible on interior of L1-ball
            inner_done = (valid_affine_step == 0);
            take_interior = ~(valid_affine_step == 1);
            lockout_dims = [];
            %things are weird on the L1 ball
        else
            [c_out, valid_affine_step] = decision_simplex(c, c_full, tau, 0); 
            lockout = (valid_affine_step == 0);
            inner_done = (valid_affine_step == 1);
            take_interior = ~(valid_affine_step == 2);
            
            %new dimensions always have a 'positive' sign, and lockout occurs
            %on a new dimension when c_aff is negative. 
            lockout_dims = find(sign(c_full) < 0 & sign(c) == 0);
        end
        
        c_aff = c_full;
        c = c_out;
        s = sign(c);
    else
        %exterior affine steps and bashing
        %must start or land on exterior
        
        %find x_aff, affine optimal point on exterior of L1-ball/simplex
        [B, p] = affine_ext_representation(s, 1, 1);
        rhs_aff = B'*(rhs_full/tau - K_full*p);
        
        if USE_CHOL
            c_aff = affine_chol_solve(B, p, K_full, L_int, rhs_full);
        else
            K_ext = B'*K_full*B; %there is probably a more efficient way to do this.
            y_aff = K_ext \ rhs_aff;
            c_aff = tau*(B*y_aff + p);
        end
        
        %clean up numerical error
        %c_aff(abs(c_aff) < 1e-9) = 0;
        
        if CS
            [c_out, valid_affine_step] = decision_l1(s, c, c_aff, tau, 1);
        else
            [c_out, valid_affine_step] = decision_simplex(c, c_aff, tau, 1); 
        end
        %an exterior step will always land on exterior. Normalize to
        %prevent propagation of floating point errors
        c_out = c_out / norm(c_out, 1) * tau; %normalize in case of numerical errors
            
        
        lockout = (valid_affine_step == 0);
        
        %new dimensions always have a 'positive' sign, and lockout occurs
        %on a new dimension when c_aff is negative. 
        if lockout
            lockout_dims = find(sign(c_aff) < 0 & sign(c) == 0);
        else
            lockout_dims = [];
        end
        
        %stable-suboptimal face state transition
        %if the face is stable, but the -gradient points inwards. Follow
        %with an interior step to leap onto a new face.
        c = c_out;
        s = sign(c);
        if valid_affine_step == 1
            grad_c = K_full * c - rhs_full;
            sg = sign(grad_c);
            if all(s(s~=0 & sg~=0) == sg(s~=0 & sg~=0))
                %stable-subopotimal face
                take_interior = 1;
            else
                %stable-optimal face
                inner_done = 1;            
            end
        end
        
        %how to define stable-suboptimal face condition?
        %ignore for now, fix later
        %TODO!
    end
    
    %slide!
    %save x_aff calls by sliding in direction of higher dimensional xaff
    slide_done = 0;
    slide_num = 0;
    if USE_SLIDE && ~inner_done && (take_interior == 0) && ~lockout
        %find a more elegant way to keep track of boundary state
        on_boundary_slide = abs(norm(c, 1) - tau) < 1e-9;
        while ~slide_done
            %find the gradient of point x
            grad_int = K_full*c - rhs_full;

            %The direction of travel is heading towards c_aff. Going directly
            %towards c_aff is impossible, so slide across the current face as
            %best as possible.
            %c_aff is a minimizer across the extended face, so by strong
            %convexity travelling towards c_aff is a valid direction     
            w_desired = c_aff - c;
            
            if CS || on_boundary_slide
                %on L1-ball or main face of simplex, use outward normals
                w_perp = (w_desired'*s)/(s'*s) * s;
                w = w_desired - w_perp;
                w(s == 0) = 0; %killed dimensions are in the nullspace of the projection
            else
                %on sidewalls, the normals are different. account for that
                %sw = zeros(size(s));
                %sw = -(c == 0); %sidewall is when c is zero
                %sw(sign(c_aff)
                %w_perp = (w_desired'*sw)/(sw'*sw) * sw;
                %w = w_desired - w_perp;
                
                %projection onto sidewalls are much simpler
                w = w_desired;
                w(c == 0) = 0;
                
            end
            %line minimization rule to find the coefficients with minimal
            %error across the line
            alpha_top = -grad_int'*w;
            alpha_bottom = w' * (K_full * w);
            alpha_slide = alpha_top/alpha_bottom;

            if isnan(alpha_slide) || (abs(alpha_bottom) < 1e-6)
                break
            end
            
            %coefficients given the line minimization rule
            c_slide = c + alpha_slide*w;
            %c_slide(abs(c_slide) <= 1e-9) = 0;
            c_slide(abs(c_slide) <= 1e-8) = 0; 
            
            if CS
                [c_slide_out, valid_slide_step] = decision_l1(s, c, c_slide, tau, 1);
            else
                [c_slide_out, valid_slide_step] = decision_simplex(c, c_slide, tau, 1); 
            end

            slide_done = (valid_slide_step == 1);
            c = c_slide_out;
            s = sign(c);
            
            %stop at a vertex, as vertices are stable
            if (nnz(c) == 1) || slide_num > 50
                slide_done = 1;
            end
            
            slide_num = slide_num + 1;
        
        end        
    end
    
    first_iter = 0;
end

c_bash = c;
on_boundary_out = abs(norm(c_bash, 1) - tau) < 1e-9; 