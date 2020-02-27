classdef chol_manager
    %CHOL_MANAGER keeps track of the cholesky decomposition of the kernel
    %matrix K. Operations include deleting columns of the kernel
    %sign-switching, and 'subtraction' going from interior to
    %exterior.
    
    
    %% Basic properties and constructor
    properties
        %original properties
        K0;         %original input kernel matrix
        rhs0;       %original right hand side
        tau;        %atomic norm constraint
        CS;         %central symmetry of atomic domain
        
        %derived properties
        L;          %current cholesky factorization of the kernel
        involved;   %indices which are nonzero in the system
        update_ext; %rank-2 update to go from exterior to interior state
    end
    
    methods
        function obj = chol_manager(K, rhs, tau, CS)
            %constructor
            %CHOL_MANAGER gets called in bash_en, before the bash_en_inner
            %loop to find stable faces. The input is the iteratively built
            %kernel and the right hand side expression.
            %
            %Input:
            %   K:      Full kernel matrix describing relations among
            %           atoms/matrix
            %   rhs:    right hand side to solve system
            %   tau:    atomic norm constraint
            %   CS:     central symmetry of system
            %
            %Output:
            %   obj:    chol_manager object (this is a constructor)
            
            %will add iterative cholesky later
            
            if nargin < 4
                CS = 0;
            end
            
            %original properties
            obj.K0 = K;
            obj.rhs0 = rhs;
            obj.tau = tau;
            obj.CS = CS;
            
            %derived properties
            N = length(rhs);
            obj.L = chol(K);
            obj.involved = (1:N)';
            obj.update_ext = [];
        end
        
        %main methods 
        
        %% conversion methods
        function obj = full_to_ext(obj, s)
            %FULL_TO_EXT convert kernels and business from interior to
            %exterior. Get it done. Input (s) is the set of signs.
            N = size(obj.L, 2);
            if nargin < 2
                s = ones(N, 1);
            end
            
            if isempty(obj.update_ext)
                obj.update_ext = obj.formulate_ext_update(obj.L, obj.involved, s);            

                %will need to replace everything with in-place stuff
                [L_out, involved_out] = obj.apply_ext_update(obj.L, obj.update_ext, obj.involved);
                obj.L = L_out;
                obj.involved = involved_out;
            else
                ME = MException('MATLAB:chol_manager:improper_conversion', ...
                                'full_to_ext called but already on exterior');                              
                throw(ME)
            end
        end
        
        function obj = ext_to_full(obj)
            %EXT_TO_FULL converts kernels and business from exterior back
            %to interior. Opposite of full_to_ext
            if isempty(obj.update_ext)
                ME = MException('MATLAB:chol_manager:improper_conversion', ...
                                'ext_to_full called but already on interior');                              
                throw(ME)
            else
                [L_out, involved_out] = obj.invert_ext_update(obj.L, obj.update_ext, obj.involved);
                obj.L = L_out;
                obj.involved = involved_out;
                obj.update_ext = [];
            end
                
        end
        
        %helper methods
        
        %% rank-1 updates and accessories
        %   adding and sign-switching indices
        
        function obj = delete_index(obj, ind)
            %DELETE_INDEX deletion of a single dimension, indexed by
            %homogenous coordinates (y, not involved/anchor)
            
            %delete index from cholesky decomposition
            %N = size(L, 2);
            obj.L = choldelete(obj.L, ind);
            
            %now kill it from all accessory data
            obj.involved(ind) = [];
            
            if ~isempty(obj.update_ext)                
                %delete the old index of the rank-2 update vectors.
                obj.update_ext.vec_pos(ind) = [];
                obj.update_ext.vec_neg(ind) = [];
                obj.update_ext.col(ind) = [];
                obj.update_ext.s(ind) = [];
            end
            
        end
        
        function obj = delete_anchor(obj)
            %DELETE_ANCHOR handles deletion of the anchor index when on
            %exterior, and remaining on exterior. Is a bit tricky to pull
            %off, need to convert to full and then go back to exterior on a
            %new anchor.
            
            %return to full space
            if isempty(obj.update_ext)
                ME = MException('MATLAB:chol_manager:improper_conversion', ...
                                'delete_anchor called but on interior');                              
                throw(ME)
            else
                s = obj.update_ext.s;
                obj = obj.ext_to_full();
                N = size(obj.L, 2);

                obj = obj.delete_index(N);

                %back onto exterior with a new anchor
                obj = obj.full_to_ext(s);            
            end
        end
        
        function obj = delete_indices(obj, drop_ind)
            %DELETE_INDICES when you want to delete multiple dimensions at
            %a time. Passing an index of '0' is dropping the anchor, and
            %this drop occurs last.
            drop_ind_s = sort(drop_ind, 'descend');
            
            for i = 1:length(drop_ind_s)
                di = drop_ind_s(i);
                
                if di == 0
                    obj = obj.delete_anchor();
                else
                    obj = obj.delete_index(di);
                end                   
            end
            
        end
        
        function L_flip = sign_switch(obj, L, s_new)
            %SIGN_SWITCH swithces the signs of atoms/kernels/pretty much
            %everything. replaces the signs of atoms a in S with -a.
            %Only applicable on a centrally symmetric atomic sets.
            %
            %Input:
            %   L:      Input cholesky decomposition
            %   s_new:  vector of signs (+- 1) describing the signs of the
            %           atoms. Just like the face of the L1-ball
            %Output:
            %   L_flip: Output cholesky decomposition. Probably should make
            %           in place
            
            %this explicit code for sign-switching is wrong. Find out why.
%             ind_flip = find(s_new < 0);            
%             L_flip = L;
%             N = size(L, 2);
%             
%             for i = 1:length(ind_flip)
%                 i_flip = ind_flip(i);
%                 
%                 %perform the flips
%                 L(1:(i_flip-1), i_flip)   = -L(1:(i_flip-1), i_flip);                
%                 L(i_flip, (i_flip+1):end) = -L(i_flip, (i_flip+1):end);
%                 
%             end


            s = ones(size(s_new));
            s(s_new == -1) = -1;
            %L_flip = L * diag(s);
            L_flip = full(bsxfun(@times,L,s'));

        end
        
                        
        %% rank-2 updates
        %   'subtraction' operation, going from interior to exterior
        function [update_ext] = formulate_ext_update(obj, L, involved, s)
            %FORMULATE_EXT_UPDATE creates rank-2 update for exterior
            %transformation of the cholesky decompsoition from full space 
            %to exterior space. Accomplishes indefinite rank-2 update by
            %two rank-one updates (one update and one downdate). Will
            %eventually replace with a rank-2 update
            %
            %Input:
            %   L:          cholesky decomposition of current kernel matrix 
            %               K (may be reduced in index by now)
            %   involved:   indices which are currently nonzero in full K
            %   s:          Signs (face of L1-ball)
            %   
            %Output:
            %   update_ext: Struct describing essential information of the
            %               rank-2 update taking the full kernel to the
            %               exterior kernel.
            N_involved = size(L, 2);
            if nargin < 3
                involved = (1:N_involved)';
            end
            
            if nargin < 4
                s = ones(N_involved, 1);
            end
            
            %sign switch if necessary
            L_s = obj.sign_switch(L, s);
            
            %find the anchor index 
            anchor = involved(end);
            
            %find the last column of the kernel matrix. This is the column
            %of the kernel associated with the anchor index, which will be
            %dropped. 
            L_col = L_s(:, end);
            K_col = L_s'*L_col;
            
            N_involved = length(K_col) - 1;
            
            %form the goldfarb update, which is actually rank-deficient
            %K_ext = (I + v*u')*K_full*(I + u*v');
            %u = zeros(N, 1);
            %u(end) = 1;
            v = -ones(N_involved, 1);
            y = K_col(1:end-1);
            
            %rho = u'*K*u = u'*y is K(anchor, anchor)
            rho = K_col(end);
%             Lu = L_s*u;
%             y = L_s'*Lu; %(y = Ku)
        
            %K_ext = K_full + v*y' + y*v' + (u'Ku)v*v'
            %where u'Ku = u'*y
            
            %K_ext = K_full + Z*B*Z'
            Z = [v y];
            B = [rho 1; 1 0];
            
            %eigendecomposition
            [V, D] = eig(B);
            
            %eigenvalues
            lam = diag(D); %one eigenvalue is negative, the other positive
            i_pos = find(lam > 0);
            i_neg = 3-i_pos;
            
            lam_pos = lam(i_pos);
            lam_neg = lam(i_neg);
            
            %eigenvectors
            if N_involved > 0
                skew_vec = Z * V;
                
                vec_pos = skew_vec(:, i_pos);
                vec_neg = skew_vec(:, i_neg);
            else
                vec_pos = [];
                vec_neg = [];
            end
            %formulate the update
            update_ext = struct;
            update_ext.anchor = anchor;
            update_ext.col = K_col;
            
            %indefinite update:
            %positive update
            update_ext.vec_pos = vec_pos;
            update_ext.lam_pos = lam_pos;
            
            update_ext.vec_neg = vec_neg;
            update_ext.lam_neg = lam_neg;
            
            update_ext.V = V;
            update_ext.D = D;
            
            
            %signs which got flipped
            update_ext.s = s(1:end-1);
            update_ext.s_anchor = s(end);
        end
        
        function [L, involved_out] = apply_ext_update(obj, L, update_ext, involved)
            %APPLY_EXT_UPDATE turns the interior L into an exterior L
            %uses the previously found update_ext from formulate_ext_update
            %
            %Input:
            %   L:          Cholesky decomposition of kernel (interior)
            %   involved:   Active indices of full (all) dimensions
            %   update_ext: Rank-2 update which changes the kernel
            %
            %Output:
            %   L_out:          Cholesky decomposition of kernel (exterior)
            %   invovled_out:   Active indices of exterior dimensions.
            %                   Excludes the anchor index (in update_ext)
            N = size(L, 2);
            if nargin < 4
                involved = (1:N)';
            end
            
            %vectors for the rank one update
            v_pos = sqrt(update_ext.lam_pos)*update_ext.vec_pos;
            v_neg = sqrt(-update_ext.lam_neg)*update_ext.vec_neg;
            
            %get rid of the anchor, which is always at the last column
            involved_out = involved(1:end-1);
            
            %drop the anchor column from L
            L = obj.sign_switch(L, [update_ext.s; update_ext.s_anchor]);
            L = choldelete(L, N);
            
            %now perform the low rank updates, which are full rank
            %(finally). Start with the update to boost the eigenvalues, and
            %then perform the downdate.
            if N > 1
                L = cholupdate(L, v_pos, '+'); 
                L = cholupdate(L, v_neg, '-');
            end
                                   
        end
        
        function [L, involved_out] = invert_ext_update(obj, L, update_ext, involved)
            %INVERT_EXT_UPDATE transitions from an exterior kernel to a
            %full kernel, inverting apply_ext_update
            
            %add the anchor index to the set of active indices
            N = size(L, 2);
            if nargin < 4
                involved = (1:N)';
            end
            
            involved_out = [involved; update_ext.anchor];
            
            %invert the rank-2 update
            v_pos = sqrt(update_ext.lam_pos)*update_ext.vec_pos;
            v_neg = sqrt(-update_ext.lam_neg)*update_ext.vec_neg;
                                    
            %Perform the low rank updates, which are full rank
            %(finally). Start with the update to boost the eigenvalues, and
            %then perform the downdate.
            if ~isempty(v_neg)
                %may have case where the anchor is the last dim standing.
                %account for this by only performing chol updates if L~=[].
                L = cholupdate(L, v_neg, '+');
                L = cholupdate(L, v_pos, '-');
                
                %add back the anchor column to L
                K_col = update_ext.col;
                K_col_side = K_col(1:end-1);
                K_col_diag = K_col(end);
                L = chol_append(L, K_col_side, K_col_diag);
            else
                %L is empty, so the only dimension remaining is the anchor.
                %K is then the number stored in update_ext.col.
                L = sqrt(update_ext.col);
            end
            
            

            L = obj.sign_switch(L, [update_ext.s; update_ext.s_anchor]);
            %strange numerical behaviour here. Until I can get a closed
            %form solution to sign switching where the diagonal stays
            %positive, the chol_append step flips signs and keeps the last
            %element positive. sign switching wil then ruin this, where the
            %element at the bottom corner stays negative. Invert this
            %afterwards, only if the sign of the anchor is negative.
            if update_ext.s_anchor < 0
                L(end, end) = -L(end, end);
            end
        end
        
        
        %% Solving system
        function r = system_rhs(obj)
            %finds the right hand side of the equations, for use in
            %solve_system. Equivalent of r = B'*(r0/tau - Kp) or r =
            %r0(involved)
            
            if isempty(obj.update_ext)
                %system on interior
                r = obj.rhs0(obj.involved)/obj.tau;
            else
                %system on exterior           
                %
                %These equations are wrong. Fix them.
                
%                 s = obj.update_ext.s;
%                 si = s(1:end-1); %sign of involved components
%                 sa = s(end);     %sign of anchor
%                 r_vec   = si.*obj.rhs0(obj.involved)/obj.tau - sa*obj.K0(obj.involved, anchor);
%                 r_const = sa*(obj.rhs0(anchor)/obj.tau - obj.K0(anchor, anchor));
%                 
%                 r = r_vec - r_const;
                
                %fix later, check correctness first
                 %sc = zeros(size(obj.rhs0));
                 

                 active = [obj.involved; obj.update_ext.anchor];

                 
                 [B, ~] = affine_ext_representation([obj.update_ext.s; obj.update_ext.s_anchor], 1, 1);
                 r = B' * (obj.rhs0(active) / obj.tau - obj.update_ext.col);
                
                 
                
%                 %hopefully fixed section
%                 s = obj.update_ext.s;
%                 active = [obj.involved; obj.update_ext.anchor];
%                 
%                 r_inner = obj.rhs0(active)/obj.tau - s(end) * obj.K0(active, obj.update_ext.anchor);
%                 r_s = s .* r_inner;
%                 r = r_s(1:end-1) - r_s(end);
%                 
                %r = (obj.rhs0(obj.involved)/obj.tau - obj.K0(obj.involved, anchor)) - (obj.rhs0(anchor)/obj.tau - obj.K0(anchor, anchor));
            end
        end
        
        function y_out = solve_system(obj)
            %SOLVE_SYSTEM find the homogenous coordinates y in the system
            %K_ext y = rhs_ext. K_ext = B'KB, where K=K0. This is the key
            %to affine bashing. 
            %
            %Input:
            %   obj:    Everything is handled by the cholesky manager.
            %           Used quantities are: L (cholesky), tau (atomic norm
            %           penalty), involved (active indices), anchor
            %           (optional).
            %
            %Output;
            %   y_out:  Homogenous coordinates describing output
            
            %Figure out the right hand side first
            %is a shortcut on the expression r = B'(r0/tau - K0 p)            
            r = obj.system_rhs();
            y_out = chol_solve(obj.L, r);
%             if isempty(obj.update_ext)
%                 %system on interior
%                 %r = obj.rhs0(obj.involved)/obj.tau;
%                 r = obj.rhs0(obj.involved);
%                 y_out = chol_solve(obj.L, r);   
%             else
%                 anchor = obj.update_ext.anchor;
%                 %system on exterior                
%                 r = (obj.rhs0(obj.involved)/obj.tau - obj.K0(obj.involved, anchor)) - (obj.rhs0(anchor)/obj.tau - obj.K0(anchor, anchor));
%                 y_out = chol_solve(obj.L, r);   
%             end
% 
%             %solve the cholesky system
                     
        end        
        
        function c = coeff(obj, y)
            %COEFF retrieve coefficients loading over the atoms. 
            %c = By + p
            %
            %Input:
            %   y:  Homogenous coordinates describing input
            %
            %Output:
            %   c:  Coefficients on the atoms
            
            c = zeros(length(obj.rhs0), 1);
            
            if isempty(obj.update_ext)
                %on interior, c = By
                By_val = obj.tau .* y;
                c(obj.involved) = By_val;
            else
                %on exterior, c = By + p
                s = obj.update_ext.s;
                By_val_involved = obj.tau * s .* y;
                c(obj.involved) = By_val_involved;
                
                s_anchor = obj.update_ext.s_anchor;
                By_val_anchor   = obj.tau * s_anchor * (1 - sum(y));
                c(obj.update_ext.anchor) = By_val_anchor;
            end
        end
        
        function y = homog(obj, c)
            %HOMOG convert from coefficients c into homogenous coordinates:
            %c = By + p.
            
%            i_active = find(c);
%             involved = i_active(1:end-1);
%             anchor   = i_active(end);
            
            %ignore signs for a minute
            c_active = c(c~=0);
            alpha_weights = sign(c_active) .* c_active / obj.tau;
            
            if isempty(obj.update_ext)
                %on interior
                y = alpha_weights;
            else
                %on exterior
                y = alpha_weights(1:end-1);
            end
            %y = alpha_weights(1:end-1) - alpha_weights(end);
        end
        
        %% Intersection code
        
        function [y_out, change, drop_ind, lockout] = intersect_int(obj, y_old, y_aff)
            %INTERSECT_INT intersection on the interior of the ball, path
            %from c_old -> c_aff in coefficients c.
            %
            %Input:
            %   y_old:      Previous value of coefficients
            %   y_aff:      New value of coefficients, minimizer across its
            %               affine space
            %
            %Output:
            %   y_out:      Feasible point as close to y_old as possible
            %   change:     Stay on interior (0) or switch to exterior (1)
            %   drop_ind:   Dimensions to drop. Only applicable on the
            %               simplex, when interior steps can lock out.
            if obj.CS
                [y_out, ~, alpha_max] = intersect_l1(y_old, y_aff, 1);
                drop_ind = [];
                lockout = 0;
            else
                [ y_out, ~, alpha_max ] = intersect_simplex(y_old, y_aff, 1);
                
                %lockout can occur on the simplex in terms of choice
                lockout_dims = find((y_old == 0) & (y_aff < 0));
                
                %sidewall intersection can also happen. Keep that in mind.
                drop_dims = find((y_old ~= 0) & (y_out == 0));
                drop_ind = [lockout_dims; drop_dims];
                lockout = (alpha_max == 0);
            end
            
            %change = (alpha_max == 1);
            %y_out always respects the constraint.
            change = (abs(norm(y_out, 1) - 1) < 1e-12);
        end
        
        function [y_out, stable_face, drop_ind, lockout] = intersect_ext(obj, y_old, y_aff)        
            %INTERSECT_EXT intersection on the exterior of the ball
            %(simplex), path from c_old -> c_aff in coefficients c.
            %
            %Input:
            %   y_old:  Previous value of coefficients
            %   y_aff:  New value of coefficients, minimizer across its
            %           affine space
            %
            %Output:
            %   y_out:          Feasible point as close to c_old as possible
            %   stable_face:    Whether the face s is a stable face (y_aff
            %                   is a feasible point)
            %   drop_ind:       If the face is unstable, drop_ind tells
            %                   what index to drop. If drop_ind = 0, drop
            %                   the anchor. 
            
            %in homogenous coordinates, everything sums to 1
            [ y_out, ~, alpha_max ] = intersect_simplex(y_old, y_aff, 1);                                            
            
            if alpha_max == 1                
                %stable face (good)
                stable_face = 1;
                lockout = 0;
                drop_ind = [];  
            else
                %unstable face
                stable_face = 0;
                
                if alpha_max == 0
                    %lockout occurs, can't move
                    %should add support for arbitrarily choosing anchor
                    lockout = 1;
                    lockout_dims = find((y_old == 0) & (y_aff < 0));
                    if isempty(lockout_dims)
                        %anchor locked out
                        lockout_dims = 0;
                    end
                    drop_ind = lockout_dims;
                else
                    %no lockout, free to move and select
                    lockout = 0;
                    lockout_dims = [];
                    
                    %sidewall intersection
                    drop_ind = find(abs(y_out) < 1e-12);
                    
                    %drop the anchor
                    if abs(sum(y_out)-1) < 1e-11
                        drop_ind = [drop_ind; 0];
                    end
                end
                %sidewall intersection
                
                

                
               % drop_ind = [lockout_dims; drop_ind_list];
                

            end
        end
        
        function [obj, y_out, stable_face, lockout, drop_ind] = single_bash(obj, y_old)
            %SINGLE_BASH single iteration of affine bashing, until a stable
            %face gets hit. Given homogenous coordinates in y_old.
            
            y_bash = obj.solve_system();
            if isempty(obj.update_ext)
                %interior
                [y_int, change, drop_ind, lockout] = obj.intersect_int(y_old, y_bash);
                
                %drop indices
                y_int(drop_ind) = [];
                obj = obj.delete_indices(drop_ind);
                
                %possibly hop to exterior and return
                if change && (length(y_int)>1)
                    obj = obj.full_to_ext(sign(y_int));
                    y_out = abs(y_int(1:end-1));
                    stable_face = 0;
                else
                    y_out = y_int;
                    stable_face = isempty(drop_ind);
                end
                    
            else
                %exterior
                [y_int, stable_face, drop_ind, lockout] = obj.intersect_ext(y_old, y_bash);
                
                %delete indices (lockout/drop dimensions/change anchors)
                if ~isempty(drop_ind) && (drop_ind(end) == 0)
                    y_int(drop_ind(1:end-1)) = []; %other dimensions
                    y_int(end) = []; %the new anchor
                else
                    y_int(drop_ind) = [];
                end
                    
                obj = obj.delete_indices(drop_ind);
                
                y_out = y_int;
            end
        end
        
        function grad = gradient_full(obj, y)
            %GRADIENT_FULL computes the interior (full) gradient at point y.
            %If the system is on the interior, this is simply Kx - rhs,
            %among the involved (nonzero) coordinates. Things are a bit
            %harder and more invovled when on the exterior. Want to be able
            %to detect stable-suboptimality without performing the rank-2
            %update to get back on interior.
            
            if isempty(obj.update_ext)
                %on interior
                Ly = obj.L * y;
                Ky = obj.L' * Ly;
                
                grad = Ky - obj.rhs0(obj.involved);
            else
                %on exterior
                
                %find the interior kernel among the involved dimensions,
                %multiplied by the signs times coordinates y
                %yf = [y; (1-sum(y))];
                s = [obj.update_ext.s; obj.update_ext.s_anchor];
                y_anc = 1 - sum(y);
                
                sy = obj.update_ext.s .* y;
                sy_anc = obj.update_ext.s_anchor * (1 - sum(y));
                
                %K_int = K_ext - lam+ v+ v+' - lam- v- v-'
                %there is a subtraction on lam- because lam- is positive.
                
                %from the exterior kernel
                Ly_ext = obj.L * y;
                Ky_ext = obj.L' * Ly_ext;
                
                %rank-2 update
                lp = obj.update_ext.lam_pos;
                vp = obj.update_ext.vec_pos;
                ln = obj.update_ext.lam_neg;
                vn = obj.update_ext.vec_neg;
                col = obj.update_ext.col;
                
                pos_ext = lp * vp * (vp' * y);
                neg_ext = ln * vn * (vn' * y);
                
                beta_involved = [(Ky_ext - pos_ext - neg_ext); 0];
                
                %now for the section determined by the anchor column
                beta_anchor = [col(1:end-1) * y_anc; (sum(col(1:end-1).*y) + col(end)*y_anc)];
                
                grad = obj.tau*s.*(beta_involved + beta_anchor) - obj.rhs0([obj.involved; obj.update_ext.anchor]);
            end
        end
        
        
    end

end
