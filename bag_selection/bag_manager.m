classdef bag_manager
    %BAG_MANAGER Class to handle bagging of 1d atoms
    %   is a property of the bash_manager
    
    properties
        %Atomic domain and properties
        norm_type;      %what norm to use
        sparse_vec;     %sparse vectors, like L1 and simplex, don't need fancy SVD low-rank updates
        tau;
        
        %Thin SVD of matrix of atoms
        U;
        Sig;
        V;
        
        %additional information (weights)
        w;        
        %only add one at a time
        FCFW;
        %group lasso + partition and magic
        small_groups;
        
        %cluster vs. subspace
        use_cluster;
        
    end
    
    methods
        function obj = bag_manager(tau, norm_type, w, FCFW)
            %BAG_MANAGER Construct an instance of this class
            %   Start out with nothing
            obj.tau = tau;
            obj.norm_type = norm_type;                       
            
            %sparse vectors are simplex, l1 ball
            obj.sparse_vec = (ischar(norm_type) && strcmp(norm_type, 'simp')) || ((isnumeric(norm_type) && norm_type == 1) && ~isstruct(w));
            %Fully Corrective Frank Wolfe, only add one atom at a time
            obj.FCFW = FCFW;
            
            
            %for testing purposes, if input is a group, then use cluster
            %bagging
            obj.small_groups = isstruct(w);
            
            obj.use_cluster = 1;
            %obj.use_cluster = obj.sparse_vec || obj.FCFW || obj.small_groups || strcmp(norm_type, 'chain');
            %for FCFW and sparse_vec, thin-svds are not needed.
            
            %SVD is empty
            obj.U = [];
            obj.Sig = [];
            obj.V = [];
            
            %weights
            if nargin < 3
                obj.w = 1;
            else
                obj.w = w;
            end
            
            
            
        end
        
        function [obj, S_bag, DG] = bag_atoms(obj, grad, x, N_bag, S0)
            %BAG_ATOMS add N_bag atoms to the system
            %Input:
            %   grad:   Gradient at x
            %   x:      Current coordinate x
            %   N_bag:  How many atoms to add if possible
            %   S0:     Original set of atoms (get rid of this)
            %   on_boundary:    If x is currently on the boundary
            %
            %Output:
            %   S_bag:  Atoms in the bag to add to system
            %   DG:     Duality gaps of each atom in the bag
            %fix this later
            if obj.use_cluster
                if obj.FCFW
                    N_bag = 1;

                    atom = LMO_1d(-grad, obj.norm_type, obj.w);                    
                    w_dir = obj.tau*atom - x;
                    DG = -grad'*w_dir;
                    if DG < 1e-4
                        S_bag = [];
                    else
                        S_bag = atom;
                    end
                else                                
                    [S_bag, DG] = bag_1d_cluster(grad, x, N_bag, obj.norm_type, obj.w, obj.tau);
                end
            else
                [S_bag, ~, ~, obj.U, obj.Sig, obj.V] = bag_1d_svd(grad, x,...
                    N_bag, S0, obj.tau, obj.norm_type, obj.w, 0, obj.U, obj.Sig, obj.V);
            end
            
        end
        
        function obj = reset_svd(obj, S)
            %reset the svd of the atoms
            [obj.U, obj.S, obj.V] = svd(S);
        end
        
        function obj = add_col(obj, add_col, ind, force_orth)
            %ADD_COL adds a new column into the thin svd
            %
            %Input:
            %   add_col:    Column to add
            %   ind:        Index at which to add the column
            %   force_orth: Force orthogonality in the low rank update
            if ~obj.use_cluster
                d = size(obj.V, 2);

                if nargin < 4
                    force_orth = 0;
                end

                if nargin < 3
                    %add to last column
                    ind = d + 1;
                end

                %perform the column addition update

                if ~isreal(obj.U)
                    add_col = complex_fold(add_col, 1);
                end
                
                [obj.U, obj.Sig, obj.V] = add_svd_column(obj.U, obj.Sig, obj.V, ...
                    add_col, force_orth, ind);
            end
        end
        
        function obj = delete_col(obj, ind, drop_col)
            if ~obj.use_cluster
                if ~isreal(obj.U)
                    drop_col = complex_fold(drop_col, 1);
                end
                [obj.U, obj.Sig, obj.V] = delete_svd_column(obj.U, obj.Sig, obj.V, ind, drop_col);
            end
        end
        
        function obj = sign_switch(obj, S, s_new)
            %SIGN_SWITCH flip things around
            [d, r] = size(S);
            

            if ~obj.use_cluster
                ind_flip = find(s_new == -1);
                num_flip = length(ind_flip);

                if num_flip > 0
                    %form the block matrices for flipping
                    
                    
                    A_flip = -2*S(:, ind_flip);
                    i_flip = ind_flip;    
                    j_flip = 1:num_flip;
                    v_flip = ones(size(i_flip));
                    B_flip = sparse(i_flip, j_flip, v_flip, r, num_flip);

                    %hack for complex numbers
                    if ~isreal(obj.U)
                        A_flip = complex_fold(A_flip, 1);
                    end
                    
                    [obj.U, obj.Sig, obj.V] = svd_update(obj.U, obj.Sig, obj.V, A_flip, B_flip, 0);
                end
            end
        end
    end
end

