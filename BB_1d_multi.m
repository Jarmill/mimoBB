function [x_final, S_final, c_final, atom_source_final, run_log] = BB_1d_multi(A, b, cons, delta)

%BB_1d_multi Bag and Bash implementation of 1d-penalized regression. All under one
%roof. Minimizes LSQ + L2 penalty (least squares + l2 penalty).
%Includes multiple atomic constraints as penalties.
%
%Supported norms (or gauge functions) are:
%   Lp norms where p >= 1 (including L1 and Linf) [p]
%   Simplex constraint (L1 + positive) ['simp']
%   Reweighted heuristic (diagonal weights in w, ||wx||_p <= tau)
%   k-chain lasso (length w) ['chain']
%   OWL norm (weights w) ['OWL']
%   Point cloud (points in w) ['point']
%   Linear system/set of stable poles ['poles']
%
%Constraints are stored in array called 'cons'
%   norm_type:  norm descripton of current penalty
%   tau:        radius of penalty
%   w:          accessory penalty information
%   index:      which indices are involved in the penalty
%
%
%Input:
%   A:          Data matrix, m x n (records by number of predictors)
%   b:          Answer vector m x 1
%   cons:       Set of atomic constraints
%   delta:      L2 regularization parameter

%Output:
%   x_final:    x that minimizes LSQ + L2 objective
%   S_final:    Atoms used to represent x_final, active set
%   c_final:    Weights on atoms in S_final to form x_final
%   run_log:    Tracking information on run

%start setup of general parameters
if nargin < 6
    w = 1;
end
%number of variables
N = size(A, 2);
x = sparse(N, 1);

%sparsify this later
%Ax = zeros(size(A, 1), 1);
Ax = cell(1, length(cons));


on_boundary = 0;

res = A*x - b;
grad = A'*res + delta*x;
error_orig = b'*b;
error_old = error_orig;
error_gap = Inf;

terminate = 0;
k = 1;


%bash management
FCFW = 1;

%BM = bash_manager(b_bash, tau, delta, norm_type, w, FCFW);
BM = multi_bash_manager(b, cons, delta, FCFW);
y = cell(length(cons), 1);

is_complex = 0;

% if BM.bagger.sparse_vec && ~is_complex
%     bag_int_increase = 8;
%     bag_ext_increase = 4;
%     N_max = Inf;
% else
%     bag_int_increase = 3;
%     bag_ext_increase = 1;
%     N_max = 10;
% end
% 
% 
% N_bag = bag_int_increase;
N_bag = 1;

%visualization
visualize = 0;
visualize_end = 0;
visualize_delay = 0;

will_visualize = visualize || visualize_end || visualize_delay;
%first_boundary = 0;


    run_log.error_list  = [];
    run_log.atomic_norm = zeros(BM.num_cons, 0);
    run_log.cardinality = [];
    run_log.num_attempted = [];
    run_log.duality_gap = [];
    run_log.time = [];
if will_visualize
    figure;
    color_list = get(gca,'ColorOrder');    
    if is_complex
        th = linspace(0, 2*pi, 51);
        circ = [cos(th); sin(th)]; 
        ind_range = 1:length(x);
        view_complex = 3;
        aspect_complex = [4 1 1];
    end
    
end

tic;
%% main bag and bash loop
while ~terminate  
    [BM, S_bag, DG] = BM.bag_atoms(grad, x, N_bag);
    
    %if the bag is empty, then no more atoms can be added to the system
    %this is a generalized termination condition
    %if isempty(S_bag) || abs(error_gap) < 1e-7    
    if all(cellfun(@isempty, S_bag))
    %if all(DG < 1e-4)
        terminate = 1;
        %DG = 0;        
        
        %input is all zeros, or the optimal point is 0
        if k == 1 
            grad_new = grad;
        end
    else
	
		%choose system based on max DG
		%in future,  explore cyclical and random.
		[DG_max, ind_cons] = max(DG);
		
		%indices of variables currently used in optimization
		ind_curr = BM.cons{ind_cons}.index;
        
		
		S_bag_curr  = S_bag{ind_cons};
		AS_bag_curr = A(:, ind_curr)*S_bag_curr;
        y_curr = y{ind_cons};
%         if is_complex
%             S_bag = complex_unfold(S_bag, 1);
%             AS_bag = complex_unfold(AS_bag, 1);
%         end
        
        %BASH
        %Find an optimal loading over current and new atoms
        [BM.basher{ind_cons}, ynew_curr] = ...
			BM.basher{ind_cons}.bash(y_curr, S_bag_curr, AS_bag_curr);
        
		y{ind_cons} = ynew_curr;
		
        %c_new = BM.basher{ind_cons}.get_c(y{ind_cons});
        
        %on_boundary_new = BM.update_ext.exterior;
        %error_new = BM.get_error(y_new);
        
        %Find statistics of output, including properties of x_new
%         N_survived = size(S_bag, 2);
%         N_dropped = 0;
%         %update the bag size
%         N_bag_next = bag_size_update(N_survived, N_dropped, bag_int_increase, bag_ext_increase, ...
%             on_boundary, on_boundary_new, N_max);
         
        x_new_curr  = BM.basher{ind_cons}.get_x(ynew_curr);
        Ax_new_curr = BM.basher{ind_cons}.get_Ax(ynew_curr);
        
        x_new = x;
        x_new(ind_curr) = x_new_curr;
        
        Ax_new = Ax;
        Ax_new{ind_cons} = Ax_new_curr;
        %S_new_curr  = BM.basher{i}.get_S();
                
        
%         if is_complex
%             x_new = complex_fold(x_new, 1);
%             Ax_new = complex_fold(Ax_new, 1);
%             S_new = complex_fold(S_new, 1);
%         end

        Ax_new_sum = sum(cell2mat(Ax_new), 2);
        
        error_new = 0.5*norm(Ax_new_sum-b)^2 + 0.5*delta*norm(x_new)^2;
        error_gap = error_old - error_new;
        
        %grad_new = BM.grad(A, y_new);
        grad_new = A'*(Ax_new_sum-b) + delta*x_new;
    end
    
    
    %% logging of the run
    %c_new = BM.get_c(y_new);
    run_log.error_list(k)  = error_new;
    %run_log.atomic_norm(:, k) = BM.tau' .* BM.get_slant_constraints(c_new);
    %atom_norm_new = BM.get_slant_constraints(c_new);
    %run_log.atomic_norm = [run_log.atomic_norm atom_norm_new];
    %run_log.cardinality(k) = nnz(c_new);
    %run_log.num_attempted(k) = size(S_bag, 2) + length(c);
    run_log.duality_gap(k) = max(DG);   
    run_log.time(k) = toc;
    
    %% visualizations go here
    
    if visualize || (terminate && visualize_end) || (visualize_delay == k)
        %the next project: putting this into a BB_visualizer class
        clf       

    %gradient (previous)
        subplot(4, 2, 1)
        hold on
        if is_complex
            scatter3(ind_range(x ~= 0), real(grad(x ~= 0)), imag(grad(x ~= 0)), 100, '.k')
            scatter3(ind_range(x == 0), real(grad(x == 0)), imag(grad(x == 0)), 100, [21 200 225]/255.0, '.')

            for i = ind_range    
                if length(w) > 1
                    grad_curr = grad(i)/w(i);
                else
                    grad_curr = grad(i);
                end
                r = abs(grad_curr);
                if x(i) ~= 0
                    plot3([i, i], [0, real(grad_curr)], [0, imag(grad_curr)], 'k')
                    if norm_type == 1
                        plot3(ones(size(th))*i, r*circ(1, :), r*circ(2, :), ':k')
                    end
                else
                    plot3([i, i], [0, real(grad_curr)], [0, imag(grad_curr)], 'Color', [21 200 225]/255.0)
                end
            end
            xlabel('$j$', 'Interpreter', 'latex')
            ylabel('Re $\nabla f(x_\mathit{prev})_j$', 'Interpreter', 'latex')
            zlabel('Im $\nabla f(x_\mathit{prev})_j$', 'Interpreter', 'latex')
            pbaspect(aspect_complex)
            view(view_complex)
        else
            stem(find(~x), grad(x == 0), '.', 'Color', [21 200 225]/255.0)
            stem(find(x), grad(x ~= 0), '.k')

            xlabel('$j$', 'Interpreter', 'latex')           
            ylabel('$\nabla f(x_\mathit{prev})_j$', 'Interpreter', 'latex')   
        end
        title('Gradient (previous)')
        hold off

    
    
    %gradient
        subplot(4, 2, 2)
        hold on
        if is_complex
            scatter3(ind_range(x_new ~= 0), real(grad_new(x_new ~= 0)), imag(grad_new(x_new ~= 0)), 100, '.k')
            scatter3(ind_range(x_new == 0), real(grad_new(x_new == 0)), imag(grad_new(x_new == 0)), 100, [21 200 225]/255.0, '.')

            for i = ind_range    
                if length(w) > 1
                    grad_curr = grad_new(i)/w(i);
                else
                    grad_curr = grad_new(i);
                end
                r = abs(grad_curr);
                if x_new(i) ~= 0
                    plot3([i, i], [0, real(grad_curr)], [0, imag(grad_curr)], 'k')
                    if norm_type == 1
                        plot3(ones(size(th))*i, r*circ(1, :), r*circ(2, :), ':k')
                    end
                else
                    plot3([i, i], [0, real(grad_curr)], [0, imag(grad_curr)], 'Color', [21 200 225]/255.0)
                end
            end                        
            xlabel('$j$', 'Interpreter', 'latex')
            ylabel('Re $\nabla f(x)_j$', 'Interpreter', 'latex')
            zlabel('Im $\nabla f(x)_j$', 'Interpreter', 'latex')
            
            pbaspect(aspect_complex)
            %view(view_complex)
            if isnumeric(norm_type) && (norm_type == 1)
                view(-90, 0)
            end
        else
            stem(find(~x_new), grad_new(x_new == 0), '.', 'Color', [21 200 225]/255.0)
            stem(find(x_new),  grad_new(x_new ~= 0), '.k')

            xlabel('$j$', 'Interpreter', 'latex')           
            ylabel('$\nabla f(x)_j$', 'Interpreter', 'latex')   
        end
        
        if terminate
            title('Gradient (final)')
        else
            title('Gradient')
        end
        hold off

            
    %computed x = Sc
        subplot(4, 2, 3)
        hold on
%         if total_variation
%             stem(cumsum(x), '.')
%             stem(cumsum(x_new), '.')
%             xlabel('$j$', 'Interpreter', 'latex')        
%             ylabel('$x_j$', 'Interpreter', 'latex')
%         else
%             if is_complex
%                 x_active = find(x_new);
%                 for i =  1:length(x_active)
%                     ind_curr = x_active(i);
%                     x_curr = x_new(ind_curr);
%                     plot3([ind_curr, ind_curr], [0, real(x_curr)], [0, imag(x_curr)], 'color', [0.8500    0.3250    0.0980])    
%                 end
%                 scatter3(x_active, real(x_new(x_active)), imag(x_new(x_active)), 100, [0.8500    0.3250    0.0980],  '.');
%                 
%                 xlabel('$j$', 'Interpreter', 'latex')        
%                 ylabel('Re $x_j$', 'Interpreter', 'latex')
%                 zlabel('Im $x_j$', 'Interpreter', 'latex')
%                 
%                 pbaspect(aspect_complex)
%                 view(view_complex)
%             else
                stem(x, '.')
                stem(x_new, '.')
                xlabel('$j$', 'Interpreter', 'latex')        
                ylabel('$x_j$', 'Interpreter', 'latex')
        %    end
        %end
        hold off
        

        title('Regressor')
        
    %Atomic Norm
        subplot(4, 2, 4)        
        plot(0:size(run_log.atomic_norm, 2), [zeros(BM.num_cons, 1) run_log.atomic_norm]')

        hold on
        for i = 1:BM.num_cons
            plot(xlim, [BM.cons{i}.tau BM.cons{i}.tau], ':k')
        end
        hold off
        
        
        xlabel('iteration')        
        ylabel('$||x||_\mathcal{A}$', 'Interpreter', 'latex')
        title('Atomic Norm')
    
    %Basis Functions
        subplot(4, 2, 5)
        hold on        
        %color_list
%         if is_complex
%             N_atoms = size(S_new, 2);
%             for ka = 1:N_atoms
%                 c_curr = c_new(ka);
%                 i_curr = find(S_new(:, ka));
%                 color_ind = mod(ka-1, length(color_list))+1;
%                 if norm_type == 1
%                     color_curr = 'k';
%                 else
%                     color_curr = color_list(color_ind, :);
%                 end
%                 for i = 1:length(i_curr)
%                     %too many indices
%                     ic = i_curr(i);
%                     S_curr = S_new(ic, ka);
%                     plot3([ic, ic], tau*c_curr*[0, real(S_curr)], tau*c_curr*[0, imag(S_curr)], 'Color', color_curr)
%                 end
%                 scatter3(i_curr, tau*c_curr*real(S_new(i_curr, ka)), tau*c_curr*imag(S_new(i_curr, ka)), 100, color_curr, '.')
%             end
%             xlabel('$j$', 'Interpreter', 'latex')        
%             ylabel('Re $S c_i$', 'Interpreter', 'latex')
%             zlabel('Im $S c_i$', 'Interpreter', 'latex')
%             
%             pbaspect(aspect_complex)
%             view(view_complex)
%         else
%             if BM.bagger.sparse_vec
%                 basis_func = x_new;
%             else
                 basis_func = S_new*diag(c_new);
%             end
            for i = 1:size(basis_func, 2)
                basis_curr = basis_func(:, i);
                active_curr = find(basis_curr);
                stem(active_curr , basis_curr(active_curr), '.')
            end
                    
            xlabel('$j$', 'Interpreter', 'latex')        
            ylabel('$x_i = S c_i$', 'Interpreter', 'latex')
%         end
        hold off
        
        title('Basis Functions')
    
    %Coefficients
        subplot(4, 2, 6)
        hold on
        stem(c_new, '.')
        hold off
        
        ylabel('$c_i$', 'Interpreter', 'latex')
        xlabel('$i$', 'Interpreter', 'latex')
        title('Coefficients')
   
    %Error
        subplot(4, 2, 7)
        semilogy(0:length(run_log.error_list), [error_orig run_log.error_list])
        
        xlabel('iteration')        
        ylabel('$f(x)$', 'Interpreter', 'latex')
        title('Function Value (error) and Duality Gap')
        %title('Function Value (error)')
        
        yyaxis right
        if terminate
            last_atom = LMO_1d(grad_new, norm_type, w);
            last_DG = real(-grad_new' * (tau*last_atom - x_new));
            semilogy(0:(length(run_log.duality_gap)), [run_log.duality_gap last_DG])
        else
            semilogy(0:(length(run_log.duality_gap)-1), run_log.duality_gap)
        end
        ylabel('DG', 'Interpreter', 'latex')
        
    %Cardinality
        subplot(4, 2, 8)
        hold on
        plot(0:length(run_log.cardinality), [0 run_log.cardinality])
        plot(run_log.num_attempted, 'xm')
        %stem(x_new, '.')
        %if first_boundary
        %        plot([first_boundary, first_boundary], ylim, ':b')
        %end
        hold off
        
        xlabel('iteration')        
        ylabel('$||c||_0$', 'Interpreter', 'latex')
        title('Cardinality')
        
        
    %look at the output
        if visualize || (visualize_delay == k)
            keyboard;
        end
    end
    
    %% prepare for new iteration
    if ~terminate
        k = k+1;
        x = x_new;
		Ax = Ax_new;
        grad = grad_new;        
        
        %y = y_new;
        error_old = error_new;
        %on_boundary = on_boundary_new;
        %N_bag = N_bag_next;
    end
end

%final output
%package and ship out
%x_final = BM.get_x(y);
x_final = x;
S_final = cellfun(@(Bi) Bi.get_S(), BM.basher);
c_final = cellfun(@(Bi, yi) Bi.get_c(yi), BM.basher, y);
atom_source_final = BM.atom_source;

% if total_variation
%     x_final = cumsum(x_new);
%     S_final = cumsum(S_new, 1);    
% end
% 
% if is_complex
%     x_final = complex_fold(x_final, 1);
%     S_final = complex_fold(S_final, 1); 
% end