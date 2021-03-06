classdef multi_bash_manager
    %MULTI_BASH_MANAGER performs bash optimization under multiple atomic
    %constraints. This is a specialization of bash_manager. The first
    %implementation of this is geared towards LPV synthesis with multiple
    %atomic constraints, but I hope to extend this to Augmented Lagrangian
    %methods.
    
    properties
        %original properties
        %A;         %data matrix (hopefully unnecessary)
        b;          %answer vector        
        cons;       %constraint set (hard regularization)
        num_cons;   %number of constraints (atomic balls)         
        delta;      %elastic net penalty        
                
        
        %optimization of each subproblem
        basher;     %performs bash on each constraint set
    end
    
    methods
        function obj = multi_bash_manager(b, cons, delta, FCFW)
            %MULTI_BASH_MANAGER Construct an instance of this class
            %   Detailed explanation goes here
            
            %Differences between bash_manager and multi_bash_manager:
            %norm constraints stored in a cell called 'cons', each of which
            %has attributes            
            %   norm_type:  norm descripton of current penalty
            %   tau:        radius of penalty
            %   w:          accessory penalty information
            %   index:      which indices are involved in the penalty
                        
            if nargin < 4
                FCFW = 0;
            end
            
            %original properties
            obj.b = b;
            obj.cons = cons;
            obj.num_cons = length(cons);
            obj.delta = delta;                                               


            %Basher and Dancer and Prancer and Vixen
            obj.basher = cell(1, obj.num_cons);
            %obj.tau = [];
            for i = 1:obj.num_cons
                cc = obj.cons{i};
                obj.basher{i} = bash_manager(b, cc.tau, delta, cc.norm_type, cc.w, FCFW);
                %obj.tau(i) = cc.tau;                
            end
            
        end
        
        %% Adding and deleting atoms
        
        function [obj, S_bag, DG] = bag_atoms(obj, grad, x, N_bag)
            %BAG_ATOMS proposes a list of new atoms to add to the system            
            %S_bag = cell(1, obj.num_cons);            
            %S_bag = sparse(length(x), N_bag*obj.num_cons);
            
            S_bag = cell(obj.num_cons, 1);
            atom_source_bag = [];
            DG =  zeros(obj.num_cons, 1);
            
            for i = 1:obj.num_cons
                idx = obj.cons{i}.index;
                grad_curr = grad(idx);
                x_curr = x(idx);
				[obj.basher{i}, S_bag_curr, DG_curr] = obj.basher{i}.bag_atoms(grad_curr, x_curr, N_bag);                
                
                %add atoms to bag
                N_added = size(S_bag_curr, 2);             
                if N_added
                    DG(i) = DG_curr;
                    S_bag{i} = S_bag_curr;
                end
            end                        
        end
                   
        
        function  x = get_x(obj, y)   
            x = cell(length (obj.cons), 1);
            for i = 1:length(obj.cons)
                x{i} = obj.basher{i}.get_x(y{i});
            end
            
            x = cell2mat(x);
        end
                  
        function  Ax = get_Ax(obj, y)   
            Ax = cell(1, length (obj.cons));
            for i = 1:length(obj.cons)
                Ax{i} = obj.basher{i}.get_Ax(y{i});
            end
            
            Ax = sum(cell2mat(Ax), 2);
        end
        
        % %% Bashing section
        % function [obj, y] = bash(obj, y_in, ind_cons, S_bag, AS_bag)
			% %Adds the atom to the system indexed by ind_cons, and then bashes accordingly.
			% %effectively FCFW.
			% y_in_curr = y_in{ind_cons};
			% [obj.basher{ind_cons}, y_out_curr] = obj.basher{ind_cons}.bash(y_in_curr, S_bag, AS_bag);
			% y{ind_cons} = y_out_curr;
		% end
        
        
    end
end

