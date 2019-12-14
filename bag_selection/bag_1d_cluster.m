function [atom_bag_list, N_added] = bag_1d_cluster(grad, x, N_bag, norm_type, w, tau)
%BAG_1D_CLUSTER performs cluster bagging on 1d input
%is really an abstraction, since cluster bagging is not practical for most
%atomic sets. Cluster bagging only works on L1, simplex. Future support
%will be added for other gauges and norms
%
%Input:
%   grad:       gradient at current x
%   x:          Current vector x = S*c, regressor of interest
%   N_bag:      number of elements that would like to be bagged
%   norm_type:  which norm or gauge to use as the atomic set
%   w:          Set of weights on reweighted heuristic, diagonal scaling
%   tau:        Atomic constraint (not used for sparse vectors)
%
%Output:
%   atom_bag_list:  Set of atoms in the bag, to be added in BASH
%   N_added:        Number of atoms in the bag

%Figure out which norm to use. For sparse vectors, will use the specialized
%bagging routine. Otherwise, will add only one atom at a time. May add
%support for more atoms in cluster bagging later, such as for point clouds
%or pos. OWL and Linf are too big and atoms are too close, and Lp is an
%infinite atomic set.
small_groups = isstruct(w);
if ~isreal(grad) && norm_type == 1 && ~small_groups
    [atom_bag_list, N_added] = bag_complex_1d(grad, x, N_bag, tau, w);
elseif isnumeric(norm_type) && norm_type == 1 && ~small_groups
    [atom_bag_list, N_added] = bag_l1(grad, x, N_bag, w);
elseif strcmp(norm_type, 'simp') && ~small_groups
    [atom_bag_list, N_added] = bag_simplex(grad, x, N_bag, w);
elseif strcmp(norm_type, 'chain') && ~small_groups
    [atom_bag_list, N_added] = bag_chain(grad, x, N_bag, tau, w);
elseif strcmp(norm_type, 'poles')
    %set of stable poles
    %fix this later for true cluster bagging
    if isfield(w, 'systems')
        %have multiple systems each contributing poles
        %should have a way to share poles between systems? That's for the
        %future maybe.
        sys_num = length(w.systems);
        %atom_bag_list_sys = sparse(length(grad), sys_num);
        atom_cell = cell(1, sys_num);
        %DG_list = zeros(sys_num, 1);
        DG = 0;
        atom_bag_list = zeros(length(grad), 1);
        
        for i = 1:sys_num
            %active indices of system
            ind_sys = w.systems{i};
            grad_sys = grad(ind_sys);
            atom_sys = LMO_1d(grad_sys, norm_type, w);
            %atom_bag_list_sys(ind_sys, i) = atom_sys;
            atom_cell{i} = atom_sys;
            
            %now find duality gap
            %I had a bug here, was measuring only against current system
            %and not against everything
            %x_sys = x(ind_sys);
            w_dir = -x;
            w_dir(ind_sys) = w_dir(ind_sys) + tau*atom_sys;
%            DG_list(i) = -grad'*w;
            DG_only = -grad'*w_dir;
            
            %only include this system in the assemblage if it helps out on
            %its own. Not sure of how well this will work.
            if DG_only > 1e-5
                atom_bag_list(ind_sys) = atom_sys;
                DG = DG + DG_only;
            end
            %DG_list(i) = -grad_sys'*(tau*atom_sys - x(ind_sys));
        end
        if any(nnz(atom_bag_list))
            N_added = 1;
        else
            N_added = 0;
            atom_bag_list = [];
        end
        %DG = grad'*(tau*atom_bag_list - x)
        %atom_bag_list = atom_bag_list_sys(:, I_sort(DG_sort > 1e-5));
         
        
        % The big screwup that I had was thinking at the same tau-penalty
        % was shared between all of the subsystems. This is not true.
        %
        % There are multiple independent penalty functions going on here:
        %   |x_1| a_1 < tau_1
        %   |x_2| a_2 < tau_2
        % different atomic norm constraints on a partition
        % The standard FW-optimization strategy would be to combine all of
        % them together, form an atom of the concatenation. That is valid
        % in general (and will be a first implementation here), but goes
        % against the theory of hierarchical sparsity in Rajiv's
        % implementation. I will need to rewrite Bash code to allow the
        % loadings over different systems with different penalties. 
        
%         [DG_sort, I_sort] = sort(DG_list, 'descend');
%         N_possible = nnz(DG_sort > 1e-5);
%         N_added = min(N_bag, N_possible);
%         I_kept = I_sort(1:N_added);
%         
%         atom_bag_list = sparse(length(grad), N_added);
%         for i = 1:N_added
%             i_curr = I_sort(i);
%             %make this more efficient
%             atom_bag_list(w.systems{i_curr}, i) = atom_cell{i_curr};
%         end

        
        
        
        
        %will add more advanced support for selecting number of atoms
        
    else
        N_added = 1;
        atom_bag_list = LMO_1d(grad, norm_type, w);
    end
else
    %for all other norms, add one at a time. Standard FW update rule
    N_added = 1;
    atom_bag_list = LMO_1d(-grad, norm_type, w);
    
    %figure out if the bag is empty
    w_dir = tau*atom_bag_list - x;
    DG = -grad'*w_dir;
    if DG < 1e-4
        atom_bag_list = [];
    end
end