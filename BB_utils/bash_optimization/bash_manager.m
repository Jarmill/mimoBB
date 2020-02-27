classdef bash_manager
    %BASH_MANAGER keeps track of the cholesky decomposition of the kernel
    %matrix K. Operations include adding and deleting atoms sign-switching, 
    %and 'subtraction' going from interior to exterior.
    
    
    %% Basic properties and constructor
    properties
        %original properties
        %A;         %data matrix (hopefully unnecessary)
        b;          %answer vector
        tau;        %atomic norm constraint
        delta;      %L2 penalty
        CS;         %central symmetry of atomic domain
        
        %derived kernel properties
        rhs;       %current right hand side
        
        %derived atom properties
        S;          %active atoms in system
        AS;         %data matrix times active atoms        
        
        
        %derived cholesky properties
%         L;          %current cholesky factorization of the kernel
%         orig_index; %index where the bag starts. orig_index=3 means 1..3 are from previous iterations, 4+ from bag.
        update_ext; %rank-2 update to go from exterior to interior state
    
        K;                      %kernel (current, interior or exterior as required)
        atom_size_max  = 200;   %preallocate for this many atoms
        %atom_size_max  = 800;   %preallocate for this many atoms
        atom_size_incr = 2;     %when the list fills up, double the size
        atom_count = 0;
        
        surviving_ind;
        
        %bagging code
        bagger;     %handles subspace bagging
    end
    
    methods
        function obj = bash_manager(b, tau, delta, norm_type, w, FCFW)
            %constructor
            %BASH_MANAGER
            
            %bagging properties
            
            %central symmetry
            if ischar(norm_type) && (strcmp(norm_type, 'point') || strcmp(norm_type, 'simp') || strcmp(norm_type, 'pos'))
                CS = 0;
            else
                CS = 1;
            end            
            
            
            if nargin < 6
                FCFW = 0;
                
                if nargin < 5
                    w = 1;
                end
            end
            
            obj.bagger = bag_manager(tau, norm_type, w, FCFW);
            
            %original properties
            obj.b = b;
            obj.tau = tau;
            obj.delta = delta;
            obj.CS = CS;
            
            %derived kernel properties
            %obj.K0 = []; %no explicit representation of kernel needed, can
            %be generated
            obj.rhs = [];
            
            %derived cholesky properties
            %obj.L = chol(K);
%             obj.L = [];
%             obj.orig_index = 0;
            %obj.update_ext = [];
            obj.update_ext.exterior = 0;
        end
        
        %main methods
        %% addition methods
        function [obj, S_bag, DG] = bag_atoms(obj, grad, x, N_bag)
            %BAG_ATOMS proposes a list of new atoms to add to the system            
            
            [obj.bagger, S_bag, DG] = obj.bagger.bag_atoms(grad, x, N_bag, obj.S(:, 1:obj.atom_count));
        end
        
        function  obj = add_atoms(obj, S_bag, AS_bag)
            %ADD_ATOMS adds atoms S_bag and AS_bag = A*S_bag to system, and
            %updates cholesky decomposition accordingly. This is how new
            %atoms get introduced into the system.
            
            
            %if these are the first atoms to be added to the bash_manager,
            %preallocate space in K, S, and AS arrays as required
            
            if ~issparse(obj.S) && isempty(obj.S)
                %kernels, atoms, and answers
                obj.S  = sparse([], [], [], size(S_bag, 1), obj.atom_size_max);
                obj.AS = sparse([], [], [], size(AS_bag, 1), obj.atom_size_max);
                %obj.K  = sparse([], [], [], obj.atom_size_max, obj.atom_size_max);
                %obj.rhs = sparse([], [], [], obj.atom_size_max, 1);
                
                obj.K = [];
                obj.rhs = [];
                
                %masks and indexing operations
                %obj.stack = (1:obj.atom_size_max)';
                %obj.mask = logical(sparse([], [], [], obj.atom_size_max, 1));
                
                obj.atom_count = 0;
                
                %rank-2 exterior update
                %obj.update_ext.Z = sparse([], [], [], obj.atom_size_max, 2);
                %obj.update_ext.Z(:, 1) = -1; %maybe fix this later?
                obj.update_ext.anchor_ind = 0;
            end
            
            %find the addition to the kernel            
            
            %bag term
            %find the number of new atoms and the indices in which they
            %will be placed.
            N_bag = size(S_bag, 2); %number of new atoms
            %bag_ind = obj.stack(1:N_bag);
            bag_ind = obj.atom_count + (1:N_bag)';
            %obj.stack = obj.stack((N_bag+1):end);                        
            
            %activate these indices on the mask
            %old_ind = find(obj.mask); %indices where there are currently active atoms
            old_ind = 1:obj.atom_count;
            
            %useful for reusing previously dropped indices
            %ind_list = [old_ind; bag_ind];
            
            %obj.mask(bag_ind) = 1;
            
            Kbag = (AS_bag'*AS_bag) + obj.delta*(S_bag'*S_bag);
            Kbag = Kbag + 1e-10*eye(size(Kbag, 2)); %slight regularization?                       
            
            %cross term
            if ~obj.update_ext.exterior
                %interior addition of atoms
                %L_int = obj.L;
                %if isempty(obj.S)
                if isempty(old_ind)
                    Kx = [];                
                    %L_new = chol(Kbag); 
                else
                    Kx = (obj.AS(:, old_ind)'*AS_bag) + obj.delta * (obj.S(:, old_ind)'*S_bag);
                    %L_new = chol_append(L_int, Kx, Kbag); 
                end
            else
                %exterior addition of atoms, this time with masking
                
                %fetch the atoms and relevant entries of the anchor
                S_anchor = obj.update_ext.S_anc;
                AS_anchor = obj.update_ext.AS_anc;
                
                %new entries of the kernel column
                y_new = AS_bag'*AS_anchor + obj.delta * S_bag'*S_anchor;
                
                %obj.update_ext.Z(bag_ind, 2) = y_new;
                Z_new = [-1*ones(size(bag_ind)) y_new];
                obj.update_ext.Z = [obj.update_ext.Z; Z_new];
                
                %rank-2 update matrix, update is K_ext = K + Z M Z'.
                M = [obj.update_ext.rho 1; 1 0];
                
                %kernel distortions (bag)
                Kbag_pre_ext = Kbag; %for debugging purposes
                Zbag = obj.update_ext.Z(bag_ind, :)*M*obj.update_ext.Z(bag_ind, :)';
                Kbag = Kbag_pre_ext + Zbag;
                
                
                %kernel distortions (cross)
                if isempty(old_ind)
                    Kx_pre_ext = [];    
                    Kx = Kx_pre_ext;
                else
                    Kx_pre_ext = (obj.AS(:, old_ind)'*AS_bag) + obj.delta * (obj.S(:, old_ind)'*S_bag);
                    
                    Zx = obj.update_ext.Z(old_ind, :)*M*obj.update_ext.Z(bag_ind, :)';
                    Kx = Kx_pre_ext + Zx;
                end
                
            end
%             else
%                 %exterior addition of atoms
%                 %see binder
%                 L_ext = obj.L;
%                 
%                 %rank-two update vectors
%                 v_old = -ones(length(L_ext), 1);
%                 y_old = obj.update_ext.K_anc_side;
%                 
%                 v_new = -ones(N_bag, 1);
%                 y_new = AS_bag'*obj.update_ext.AS_anchor + obj.delta * S_bag'*obj.update_ext.S_anchor;
%                                 
%                 Z_new = [v_new y_new];
%                 skew_vec_new = Z_new * obj.update_ext.V;
%                 vec_pos_new = skew_vec_new(:, obj.update_ext.ind_pos);
%                 vec_neg_new = skew_vec_new(:, 3-obj.update_ext.ind_pos);
%                 
%                 Kbag_ext = Kbag + [v_new y_new] * obj.update_ext.M * [v_new'; y_new'];
%                 
%                 if isempty(obj.S)
%                     %Kx_ext = [];
%                     L_new = chol(Kbag_ext);
%                 else
%                     Kx = (obj.AS'*AS_bag) + obj.delta * (obj.S'*S_bag);
%                     Kx_ext = Kx + [v_old y_old] * obj.update_ext.M * [v_new'; y_new'];
%                     L_new = chol_append(L_ext, Kx_ext, Kbag_ext);
%                 end
%                 
%                 obj.update_ext.K_anc_side = [y_old; y_new];
%                 obj.update_ext.vec_pos = [obj.update_ext.vec_pos; vec_pos_new];
%                 obj.update_ext.vec_neg = [obj.update_ext.vec_neg; vec_neg_new];
%             end

                
            %side term
            rhs_new = AS_bag'*obj.b;
            
            %prepare output
            %obj.L  = L_new;
            obj.atom_count = obj.atom_count + N_bag;
            
            %add to atoms
            obj.S(:, bag_ind) = S_bag;
            obj.AS(:, bag_ind) = AS_bag;
            
            %add to kernel
            %cross terms
            if ~isempty(old_ind)
                obj.K(old_ind, bag_ind) = Kx;
                obj.K(bag_ind, old_ind) = Kx';
            end
            %new terms (bag)
            obj.K(bag_ind, bag_ind) = Kbag;
            
            %add to answer vector
            obj.rhs(bag_ind, 1) = rhs_new;
        end
        
        %% conversion methods
        function obj = full_to_ext(obj)
            %FULL_TO_EXT convert kernels and business from interior to
            %exterior. Get it done. Input (s) is the set of signs.
            if ~obj.update_ext.exterior
                %update_ext = obj.formulate_ext_update(obj.L);            

                %will need to replace everything with in-place stuff
                %obj = obj.apply_ext_update(update_ext);
                obj = obj.apply_ext_update_kernel();
                %obj.orig_index = obj.orig_index - 1;
            else
                ME = MException('MATLAB:chol_manager:improper_conversion', ...
                                'full_to_ext called but already on exterior');                              
                throw(ME)
            end
        end
        
        function obj = ext_to_full(obj)
            %EXT_TO_FULL converts kernels and business from exterior back
            %to interior. Opposite of full_to_ext
            if ~obj.update_ext.exterior
                ME = MException('MATLAB:chol_manager:improper_conversion', ...
                                'ext_to_full called but already on interior');                              
                throw(ME)
            else
                %obj = obj.invert_ext_update();
                obj = obj.invert_ext_update_kernel();
                %obj.orig_index = obj.orig_index + 1;
                %obj.L = L_out;
                %obj.involved = involved_out;
                %obj.update_ext = [];
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
            %obj.L = choldelete(obj.L, ind);
            
            %if ind < obj.orig_index
            %    obj.orig_index = obj.orig_index - 1;
            %end
            
            %%now kill it from all accessory data
            
            %delete the old index of atom parameters
            
            S_drop = obj.S(:, ind);
            %drop from the bagger
            obj.bagger = obj.bagger.delete_col(ind, S_drop);
            
            %drop from the basher
            %obj.S(:, ind) = [];
            %obj.AS(:, ind) = [];
            obj.rhs(ind) = [];
            
            %masking: mark that dropped index is inactive, then add that
            %position to the stack
            %obj.mask(ind) = 0;
            %obj.stack = [ind; obj.stack];
            
            %pulling: move over everything in S and AS
            obj.S(:, ind:(obj.atom_count-1)) = obj.S(:, (ind+1):obj.atom_count);
            obj.S(:, obj.atom_count) = 0;
            
            obj.AS(:, ind:(obj.atom_count-1)) = obj.AS(:, (ind+1):obj.atom_count);
            obj.AS(:, obj.atom_count) = 0;
                        
            obj.surviving_ind(ind + obj.update_ext.exterior) = [];
            
            obj.atom_count = obj.atom_count - 1;            
            
            obj.K(ind, :) = [];
            obj.K(:, ind) = []; %inefficient?
            
            if obj.update_ext.exterior
                obj.update_ext.Z(ind, :) = [];
            end
            
            %drops are cheap in the masking scheme! (without cholesky)
            
        end
        
        function obj = delete_anchor(obj)
            %DELETE_ANCHOR handles deletion of the anchor index when on
            %exterior, and remaining on exterior. Is a bit tricky to pull
            %off, need to convert to full and then go back to exterior on a
            %new anchor.
            
            %return to full space
            if ~obj.update_ext.exterior
                ME = MException('MATLAB:chol_manager:improper_conversion', ...
                                'delete_anchor called but on interior');                              
                throw(ME)
            else
                obj = obj.ext_to_full();
                %N = size(obj.L, 2);

                obj = obj.delete_index(obj.update_ext.anchor_ind);

                obj = obj.full_to_ext();            
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
                    %mask_ind_drop = mask_indices(di);
                    %obj = obj.delete_index(mask_ind_drop);                    
                    obj = obj.delete_index(di);
                end                   
            end
            
        end
        
        function obj = sign_switch(obj, s_new)
            %SIGN_SWITCH swithces the signs of atoms/kernels/pretty much
            %everything. replaces the signs of atoms a in S with -a.
            %Only applicable on a centrally symmetric atomic sets.
            %
            %Input:
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
%             end


            s = ones(size(s_new));
            s(s_new == -1) = -1;

            if any(s_new == -1)
                %can have negatives, still need to make an explicit form
                
                %switch on the bagger
                obj.bagger = obj.bagger.sign_switch(obj.S, s);
                
                %switch on the basher
                %obj.L  = full(bsxfun(@times,obj.L,s'));                
                obj.K = bsxfun(@times, s, bsxfun(@times, obj.K, s'));

                obj.S(:, 1:obj.atom_count)  = bsxfun(@times,obj.S(:, 1:obj.atom_count),s');
                obj.AS(:, 1:obj.atom_count) = bsxfun(@times,obj.AS(:, 1:obj.atom_count),s');
                obj.rhs = obj.rhs .* s;
            end

        end
        
                        
        %% rank-2 updates
        %   'subtraction' operation, going from interior to exterior
        function [update_ext] = formulate_ext_update(obj, L)
            %FORMULATE_EXT_UPDATE creates rank-2 update for exterior
            %transformation of the cholesky decompsoition from full space 
            %to exterior space. Accomplishes indefinite rank-2 update by
            %two rank-one updates (one update and one downdate). Will
            %eventually replace with a rank-2 update
            %
            %Input:
            %   L:          cholesky decomposition of the kernel
            %   ind_anc:    anchor index            
            %   
            %Output:
            %   update_ext: Struct describing essential information of the
            %               rank-2 update taking the full kernel to the
            %               exterior kernel.
            %N_involved = size(L, 2);

            %anchor index is the last
            
            %find the last column of the kernel matrix. This is the column
            %of the kernel associated with the anchor index, which will be
            %dropped.             
            
            %form the goldfarb update, which is actually rank-deficient
            %K_ext = (I + v*u')*K_full*(I + u*v');
            %u = zeros(N, 1);
            %u(end) = 1;
            v = -ones(N_involved, 1);
            y = K_col(1:end-1);
            
            %rho = u'*K*u = u'*y is K(anchor, anchor)
            rho = K_col(end);
        
            %K_ext = K_full + v*y' + y*v' + (u'Ku)v*v'
            %where u'Ku = u'*y
            
            %K_ext = K_full + Z*B*Z'
            Z = [v y];
            M = [rho 1; 1 0];
            
            %eigendecomposition
            [V, D] = eig(M);
            
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
            update_ext.K_anc_side = y;
            update_ext.K_anc_diag = rho;
            
            %indefinite update:
            %positive update
            update_ext.vec_pos = vec_pos;
            update_ext.lam_pos = lam_pos;
            
            update_ext.vec_neg = vec_neg;
            update_ext.lam_neg = lam_neg;
            
            update_ext.ind_pos = i_pos;
            
            %matrix entries
            update_ext.M = M;
            update_ext.V = V;
            update_ext.D = D;
        end
        
        function obj = apply_ext_update(obj, update_ext)
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
            N = size(obj.L, 2);

            
            %vectors for the rank one update
            v_pos = sqrt(update_ext.lam_pos)*update_ext.vec_pos;
            v_neg = sqrt(-update_ext.lam_neg)*update_ext.vec_neg;
            
            
            %drop the anchor column from L
            %L = obj.sign_switch(L, [update_ext.s; update_ext.s_anchor]);
            obj.L = choldelete(obj.L, N);
            
            %now perform the low rank updates, which are full rank
            %(finally). Start with the update to boost the eigenvalues, and
            %then perform the downdate.
            if N > 1
                obj.L = cholupdate(obj.L, v_pos, '+'); 
                obj.L = cholupdate(obj.L, v_neg, '-');
            end
            
            obj.update_ext = update_ext;
            
            %clean up atoms
            obj.update_ext.S_anchor = obj.S(:, end);            
            obj.update_ext.AS_anchor = obj.AS(:, end);
            obj.update_ext.rhs_anchor = obj.rhs(end);
            
            %drop last column from the bagger
            obj.bagger = obj.bagger.delete_col(N, obj.S(:, end));
            
            obj.S(:, end) = [];
            obj.AS(:, end) = [];
            obj.rhs(end) = [];
            
            obj.update_ext.exterior = 1;
                                   
        end
        
        function obj = apply_ext_update_kernel(obj, anchor_ind)
            %APPLY_EXT_UPDATE_MASK performs transitions the system onto the
            %exterior treating index anchor_ind as the anchor. This input
            %is optional, and by default will be the active index with
            %the lowest index (a change from non-masking).
            
            %if nargin < 2 || obj.mask(anchor_ind) ~= 1
            %    anchor_ind = find(obj.mask, 1);
            %end
            
            %always use the first atom as the anchor. makes things easy.
            anchor_ind = 1;
            
            obj.update_ext.rhs_anc = obj.rhs(anchor_ind);
            obj.update_ext.S_anc = obj.S(:, 1);            
            obj.update_ext.AS_anc = obj.AS(:, 1);
            
            obj.S(:, 1:obj.atom_count-1) = obj.S(:, 2:obj.atom_count);
            obj.S(:, obj.atom_count) = 0;
            obj.AS(:, 1:obj.atom_count-1) = obj.AS(:, 2:obj.atom_count);
            obj.AS(:, obj.atom_count) = 0;
        
            obj.rhs(anchor_ind) = [];
            
            obj.update_ext.exterior = 1;
            obj.update_ext.anchor_ind = anchor_ind;
            
            %copy the anchor's kernel column over to Z
            anc_col = obj.K(:, anchor_ind);
            obj.update_ext.rho = anc_col(anchor_ind);
            anc_col(anchor_ind) = [];
            
            obj.update_ext.Z = [-ones(size(anc_col)) anc_col];
            %obj.update_ext.Z(obj.mask, 2) = obj.K(obj.mask, anchor_ind);
            
            %also found in obj.update_ext.Z(obj.update_ext.anchor_ind, 2)

            M = [obj.update_ext.rho, 1; 1, 0];
            
            %drop the anchor from the mask, then apply the rank-2 update
            
            
            %Zx = obj.update_ext.Z(obj.mask, :);
            obj.K(:, anchor_ind) = [];
            obj.K(anchor_ind, :) = [];
            if ~isempty(obj.K)
                obj.K = obj.K + obj.update_ext.Z * M * obj.update_ext.Z';
            end
            
            %all done!
            obj.atom_count = obj.atom_count - 1;
        end
        
        function obj = invert_ext_update(obj)
            %INVERT_EXT_UPDATE transitions from an exterior kernel to a
            %full kernel, inverting apply_ext_update
            
            %add the anchor index to the set of active indices
            %N = size(obj.L, 2);
            
            
            %invert the rank-2 update
            v_pos = sqrt(obj.update_ext.lam_pos)*obj.update_ext.vec_pos;
            v_neg = sqrt(-obj.update_ext.lam_neg)*obj.update_ext.vec_neg;
                                    
            %Perform the low rank updates, which are full rank
            %(finally). Start with the update to boost the eigenvalues, and
            %then perform the downdate.
            if ~isempty(v_neg)
                %may have case where the anchor is the last dim standing.
                %account for this by only performing chol updates if L~=[].
                obj.L = cholupdate(obj.L, v_neg, '+');
                obj.L = cholupdate(obj.L, v_pos, '-');
                
                %add back the anchor column to L
                %K_col = obj.update_ext.col;
                %K_col_side = K_col(1:end-1);
                %K_col_diag = K_col(end);
                K_col_side = obj.update_ext.K_anc_side;
                K_col_diag = obj.update_ext.K_anc_diag;
                
                obj.L = chol_append(obj.L, K_col_side, K_col_diag);
                
                %append anchor to atoms and properties
                obj.S   = [obj.S    obj.update_ext.S_anchor];
                obj.AS  = [obj.AS   obj.update_ext.AS_anchor];
                obj.rhs = [obj.rhs; obj.update_ext.rhs_anchor];
            else
                %L is empty, so the only dimension remaining is the anchor.
                %K is then the number stored in update_ext.col.
                obj.L   = sqrt(obj.update_ext.K_anc_diag);
                obj.S   = obj.update_ext.S_anchor;
                obj.AS  = obj.update_ext.AS_anchor;
                obj.rhs = obj.update_ext.rhs_anchor;
            end
            
            %update the bagger
            obj.bagger = obj.bagger.add_col(obj.update_ext.S_anchor);
            
            obj.update_ext = [];
                        
            obj.update_ext.exterior = 0;
            
        end
        
        function obj = invert_ext_update_kernel(obj)
            %INVERT_EXT_UPDATE_KERNEL transitions the system off the
            %exterior back to the full (possibly signed) space.            
            
            %invert the rank-2 update
            M = [obj.update_ext.rho, 1; 1, 0];
           
            obj.S(:, 2:obj.atom_count+1) = obj.S(:, 1:obj.atom_count);
            obj.S(:, 1) = obj.update_ext.S_anc;
                       
            obj.AS(:, 2:obj.atom_count+1) = obj.AS(:, 1:obj.atom_count);
            obj.AS(:, 1) = obj.update_ext.AS_anc;
                        
            obj.rhs = [obj.update_ext.rhs_anc; obj.rhs];
            obj.atom_count = obj.atom_count + 1;
            
            
            obj.K = obj.K - obj.update_ext.Z * M * obj.update_ext.Z';
            %now add back the anchor column;
            %anc_ind = obj.update_ext.anchor_ind;
            %Z_pre = obj.update_ext.Z(1:anc_ind-1 , 2);
            %Z_post = obj.update_ext.Z(anc_ind+1:end , 2);
            Z_col = obj.update_ext.Z(:, 2);
            %rely on the fact that the anchor is always at the beginning
            obj.K = [obj.update_ext.rho Z_col'; Z_col obj.K];
            %obj.K = [[obj.K(1:anc_ind-1, 1:anc_ind-1) Z_pre obj.K(1:anc_ind-1, anc_ind+1:end)],
            %        [Z_pre'; obj.update_ext.rho; Z_post'];
            %        [obj.K(anc_ind+1:end, 1:anc_ind-1); Z_post; obj.K(anc_ind+1:end, 1:anc_ind-1)]];
            %switch off exterior flag, and add anchor index back to mask
            obj.update_ext.exterior = 0;
            
            %reverted!
            
        end
        
        
        
        %% Solving system
        function r = system_rhs(obj)
            %finds the right hand side of the equations, for use in
            %solve_system. Equivalent of r = B'*(r0/tau - Kp) or r =
            %r0(involved)
            
            if ~obj.update_ext.exterior
                %system on interior
                %r = obj.rhs0(obj.involved)/obj.tau;
                r = obj.rhs / obj.tau;
            else
                %system on exterior           
                
                %These equations are wrong. Fix them.
                
%                 s = obj.update_ext.s;
%                 si = s(1:end-1); %sign of involved components
%                 sa = s(end);     %sign of anchor
%                 r_vec   = si.*obj.rhs0(obj.involved)/obj.tau - sa*obj.K0(obj.involved, anchor);
%                 r_const = sa*(obj.rhs0(anchor)/obj.tau - obj.K0(anchor, anchor));
%                 
%                 r = r_vec - r_const;
                
                %fix later, check correctness first                 
                %active = [obj.involved; obj.update_ext.anchor];
                 
                %[B, ~] = affine_ext_representation([obj.update_ext.s; obj.update_ext.s_anchor], 1, 1);                
                %r = B' * (obj.rhs / obj.tau - obj.update_ext.col);
                %r0 = [obj.rhs; obj.update_ext.rhs_anchor] / obj.tau - obj.update_ext.col;
                %r = r0(1:end-1) - r0(end);
                
                %r_side = obj.rhs/obj.tau - obj.update_ext.K_anc_side;
                %r_diag = obj.update_ext.rhs_anchor/obj.tau - obj.update_ext.K_anc_diag;
                
                
                
                %Z is [column of -1, anchor column of kernel (interior)]
                if obj.atom_count
                    r_side = obj.rhs/obj.tau - obj.update_ext.Z(:, 2);
                else
                    r_side = 0;
                end
                r_diag = obj.update_ext.rhs_anc/obj.tau - obj.update_ext.rho;
                
                
                r = r_side - r_diag;
                
                
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
            
            if obj.atom_count
                r = obj.system_rhs();
            
                y_out = obj.K \ r;
            else
                y_out = [];
            end
            
            %y_out = chol_solve(obj.L, r);
%             if ~obj.update_ext.exterior
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
            
            if ~obj.update_ext.exterior
                %on interior, c = By
                By_val = obj.tau .* y;
                c(obj.involved) = By_val;
            else
                %on exterior, c = By + p
                %s = obj.update_ext.s;
                s = 1;
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
            
            if ~obj.update_ext.exterior
                %on interior
                y = alpha_weights;
            else
                %on exterior
                y = alpha_weights(1:end-1);
            end

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
        
        %% Gradients and coefficients (access methods)
        function grad = gradient_full(obj, y)
            %GRADIENT_FULL computes the interior (full) gradient at point y.
            %If the system is on the interior, this is simply Kx - rhs,
            %among the involved (nonzero) coordinates. Things are a bit
            %harder and more invovled when on the exterior. Want to be able
            %to detect stable-suboptimality without performing the rank-2
            %update to get back on interior.
            
            if isempty(y) && obj.atom_count == 0
                %on a vertex
                grad = [];
            else            
                %on a face
                if ~obj.update_ext.exterior
                    %on interior
                    %Ly = obj.L * y;
                    %Ky = obj.L' * Ly;
                    Ky = obj.K*y;
                    grad = Ky - obj.rhs/obj.tau;
                else
                    %on exterior

                    %find the interior kernel among the involved dimensions,
                    %multiplied by the signs times coordinates y
                    %yf = [y; (1-sum(y))];
                    y_anc = 1 - sum(y);                

                    %K_int = K_ext - lam+ v+ v+' - lam- v- v-'
                    %there is a subtraction on lam- because lam- is positive.

                    %from the exterior kernel
                    %Ly_ext = obj.L * y;
                    %Ky_ext = obj.L' * Ly_ext;
                    Ky_ext = obj.K*y;
                    
                    %rank-2 update
%                     lp = obj.update_ext.lam_pos;
%                     vp = obj.update_ext.vec_pos;
%                     ln = obj.update_ext.lam_neg;
%                     vn = obj.update_ext.vec_neg;
%                     %col = obj.update_ext.col;
%                     colx = obj.update_ext.K_anc_side;
%                     cold = obj.update_ext.K_anc_diag;
% 
%                     pos_ext = lp * vp * (vp' * y);
%                     neg_ext = ln * vn * (vn' * y);
% 
%                     beta_involved = [(Ky_ext - pos_ext - neg_ext); 0];
                    
                    %K y = K_ext y - Z M Z' y
                    %Z_curr = obj.update_ext.Z(obj.mask, 2);
                    %Zx = obj.update_ext.Z;
                    
                    colx = obj.update_ext.Z(:, 2);
                    cold = obj.update_ext.rho;
                    
                    M = [cold, 1; 1, 0];
                    
                    rank2_y = obj.update_ext.Z * (M * (obj.update_ext.Z'*y));                                        
                    
                    beta_involved = [0; (Ky_ext - rank2_y)];

                    %now for the section determined by the anchor column
                    beta_anchor = [(sum(colx.*y) + cold*y_anc); colx* y_anc];

                    grad = obj.tau*(beta_involved + beta_anchor) - [obj.update_ext.rhs_anc; obj.rhs];

                end
            end
        end
        
        function obj = set_b(obj, b_in)
            %Dangerous, but this is necessary for multi-bash management
            obj.b = b_in;
            if ~isempty(obj.AS)
                obj.rhs = obj.AS(:, 1:obj.atom_count)'*obj.b;
            end
            
            if obj.update_ext.exterior
                obj.update_ext.rhs_anc = obj.update_ext.AS_anc'*obj.b;
            end
        end
        
        %need a way to handle an output of 0.
        function c = get_c(obj, y)
            if ~obj.update_ext.exterior
                %interior
                c = y;
            else
                c = [y; 1-sum(y)];                
            end
        end
        
        function x = get_x(obj, y)
            %GET_X returns x = Sc
            
            if ~obj.update_ext.exterior
                %interior
                x = obj.tau * obj.S(:, 1:obj.atom_count) * y;
            else
                if isempty(y)
                    x = obj.tau * obj.update_ext.S_anc;
                else
                    x = obj.tau * (obj.S(:, 1:obj.atom_count) * y + obj.update_ext.S_anc * (1-sum(y)));
                end
                                
            end
        end
        
        function S = get_S(obj)
            %GET_S returns S
            %for reporting purposes, includes the anchor in its proper position.            
            S = obj.S(:, 1:obj.atom_count);
            if obj.update_ext.exterior
                S = [obj.update_ext.S_anc S];
            end
            
            
        end
        
        function Ax = get_Ax(obj, y)
            %GET_X returns Ax = ASc
            
            if ~obj.update_ext.exterior
                %interior
                Ax = obj.tau * obj.AS(:, 1:obj.atom_count) * y;
            else
                %Ax = obj.AS * y;
                if isempty(y)
                    Ax = obj.tau * obj.update_ext.AS_anc;
                else
                    inv_prod = obj.AS(:, 1:obj.atom_count) * y;
                    anc_prod = obj.update_ext.AS_anc* (1-sum(y));
                    Ax = obj.tau * (inv_prod + anc_prod);
                end
                
                %Ax = Ax + obj.AS_anchor * (1-sum(y));
            end
        end
        
        function grad_all = grad(obj, A, y)
            %GRAD find the gradient of f(x) given the data matrix A and
            %homogenous coordinates y
            %
            %or maybe it would be better to have the data matrix as a
            %property of the bash_manager? I don't know.
            x = obj.get_x(y);
            Ax = obj.get_Ax(y);
            
            grad_all = A'*(Ax-obj.b) + obj.delta*x;
        end
        
        function error = get_error(obj, y)
            %GET_ERROR finds the elastic net objective error
             if  ~obj.update_ext.exterior && isempty(y) 
                 error = 0.5*sum(obj.b.^2);
             else
                 x = obj.get_x(y);
                 Ax = obj.get_Ax(y);
                 
                 error = 0.5*norm(Ax-obj.b)^2 + 0.5*obj.delta*norm(x)^2;
             end
        end
        
        %% Bashing section
        function [obj, y_out, stable_face, lockout, drop_ind] = single_bash(obj, y_old)
            %SINGLE_BASH single iteration of affine bashing, until a stable
            %face gets hit. Given homogenous coordinates in y_old.
            
            y_bash = obj.solve_system();
            if ~obj.update_ext.exterior
                %interior
                [y_int, change, drop_ind, lockout] = obj.intersect_int(y_old, y_bash);
                
                %drop indices
                y_int(drop_ind) = [];
                obj = obj.delete_indices(drop_ind);
                
               
                
                %possibly hop to exterior and return
                %if change && (length(y_int)>1)
                if change && (length(y_int)>1)
                    %sign switch if necessary
                    obj = obj.sign_switch(sign(y_int));
                    
                    obj = obj.full_to_ext();
                    %y_out = abs(y_int(1:end-1));
                    y_out = abs(y_int);
                    y_out(obj.update_ext.anchor_ind) = [];
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
                    %y_int(end) = []; %the new anchor
                    y_int(1) = []; %the new anchor
                else
                    y_int(drop_ind) = [];
                end                                
                    
                obj = obj.delete_indices(drop_ind);
                
                %vertices are stable
                if isempty(y_int)
                    stable_face = 1;
                end
                
                y_out = y_int;
            end
        end
        
        function [obj, y_new, drop_ind, bash_iter] = bash_asqp(obj, y_old)
             %use ASQP code to perform the bash step  
             
             %set up the system            
             %r = obj.system_rhs();
             r = obj.rhs;
             new_atom_added = 1;
             param.max_iter = 500;
             param.debug_mode = 0;
             param.epsilon = 1e-6;
             %param.rho = obj.tau;
             param.rho = obj.tau;
             
              %solve the problem
             [Y_new, A, bash_iter] = asqp_constrained(obj.K, r, obj.tau*y_old, param, new_atom_added);
             
             y_new = Y_new/obj.tau;
             %delete bad atoms 
             drop_ind = find(A==0);      
             obj = obj.delete_indices(drop_ind);
        end 
        
        function [obj, y] = bash(obj, y_in, S_bag, AS_bag)
            %BASH performs the core affine-bashing loop, running
            %single_bash multiple times. Finds a stable-optimal face of the
            %simplex/L1 ball.
            %
            %Input:
            %   y_in:   previous set of homogenous coordinates, starting
            %           point of bashing loop
            %   S_bag:  New atoms from the bag
            %   AS_bag: Data matrix times new atoms from the bag
            %
            %Output:
            %   y_bash: Output homogenous coordinates
            
            
            if nargin == 4
                %add new atoms
                %N_prev = length(y_in); %also nnz(obj.mask)
                N_bag = size(S_bag, 2);
                %[obj, ind_list] = obj.add_atoms(S_bag, AS_bag);
                obj = obj.add_atoms(S_bag, AS_bag);
                
                
                %adding new entries isn't this simple anymore
                y = [y_in; zeros(N_bag, 1)];
                
                %find out where the previously added indices end up
                %[~, ind_perm] = sort(ind_list);
                %destination = find(ind_perm <= N_prev);
                
                %create an expanded vector y to put the previous weights in
                %the expanded position. Not sure if there is more efficient
                %way to do this. Don't think this is the bottleneck, also
                %the list is nearly sorted (merge), might write merge code
                %y = sparse(destination, ones(N_prev, 1), y_in, N_prev+N_bag, 1);
            else
                %keep everything the same
                y = y_in;
                if ~obj.atom_count
                    return
                end
            end
            
            %debugging, tracking the error
            %error_orig = obj.get_error(y);
            %error_old = error_orig;
            
            inner_done = 0;
            bash_iter = 1;
            
            obj.surviving_ind = (1:(obj.atom_count + obj.update_ext.exterior))';
            
            
            ASQP = 0;
            if ASQP 
                %this is wrong, ASQP has an equality constraint that sum(y)
                %= rho
                [obj, y, drop_ind, bash_iter] = bash_asqp(obj, y);
            else

                %main bashing loop
                while ~inner_done
                    [obj, y, stable_face, lockout, drop_ind] = obj.single_bash(y);
                    %error_new = obj.get_error(y);
                    %error_gap = error_old - error_new;
                    if stable_face
                        if ~obj.update_ext.exterior
                            %interior, done regardless
                            inner_done = 1;
                        else
                            %exterior, need to check stable-suboptimality
                            grad_full = obj.gradient_full(y);
                            if all(sign(grad_full) == -1)                        
                                %stable-optimal face
                                inner_done = 1;
                            else
                                %stable-suboptimal face
                                %return to interior (full set of dimensions)
                                y =  [y; 1-sum(y)];
                                obj = obj.ext_to_full();
                            end
                        end
                    end
                    bash_iter = bash_iter + 1;                
                    %error_old = error_new;
                end                   
            end
            
            %return
        end
        
        
    end

end
