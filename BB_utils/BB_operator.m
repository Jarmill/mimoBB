function [x_final, S_final, c_final, run_log] = BB_operator(A, At, b, opt)
%BB_1d Bag and Bash implementation of 1d-penalized regression. All under one
%roof. Minimizes LSQ + L2 penalty (least squares + l2 penalty).
%
%Supported norms are:
%   Lp norms where p >= 1 (including L1 and Linf) [p]
%   Simplex constraint (L1 + positive) ['simp']
%   Reweighted heuristic (diagonal weights in w, ||wx||_p <= tau)
%   k-chain lasso (length w) ['chain']
%   OWL norm (weights w) ['OWL']
%   Point cloud (points in w) ['point']
%
%Input:
%   A:          Data operator (generalizes data matrix)
%   At:         Adjoint data operator
%   b:          Answer vector m x 1
%opt fields:
%   tau:        Radius of L1/ball or height of simplex
%   delta:      L2 regularization parameter
%   norm_type:  Desired regularization norm.
%   w:          Weights for reweighted heuristic, or additional information
%   num_var:    Number of variables
%
%Output:
%   x_final:    x that minimizes LSQ + L2 objective
%   S_final:    Atoms used to represent x_final, active set
%   c_final:    Weights on atoms in S_final to form x_final
%   run_log:    Tracking information on run

%start setup of general parameters
tau = opt.tau;
delta = opt.delta;
w = opt.w;
N = opt.num_var;
norm_type = opt.norm_type;


if isnumeric(A)
    %can probably do this more elegantly
    A0 = A;
    At = @(r) A0'*r;
    A  = @(x) A0*x;
    A_matrix = 1;
else
    A_matrix = 0;
end


terminate = 0;
k = 1;


%complex lasso
is_complex = opt.is_complex;

if (ischar(norm_type) && strcmp(norm_type, 'chain') && w == 1)
    norm_type = 1;
end

%bash management
FCFW = 1;
if is_complex
    %b_bash = [real(b); imag(b)];
    b_bash = complex_unfold(b, 1);
else
    b_bash = b;
end





if isfield(opt, 'warm_start')
    BM =  opt.warm_start.bash_manager;

    y_warm  = opt.warm_start.y; 
    [BM, y] = BM.bash(y_warm);
    
    c  = BM.get_c(y);
    x  = BM.get_x(y);
    Ax = BM.get_Ax(y);
    
    on_boundary = BM.update_ext.exterior;
    
    atomic_norm_start = tau*sum(c);
    cardinality_start = length(c);    
else
    y = [];
    c = [];
    x = sparse(N, 1);
    Ax = A(x); 

    BM = bash_manager(b_bash, tau, delta, norm_type, w, FCFW);
            
    on_boundary = 0;   
    
    atomic_norm_start = 0;
    cardinality_start = 0;   
end

res = Ax - b;
grad = At(res) + delta*x;
[n, nlist] = anorm_1d(x,  opt.norm_type, opt.w);

error_orig = 0.5*norm(Ax-b)^2 + 0.5*delta*norm(x)^2;
error_old = error_orig;
error_gap = Inf;

if BM.bagger.sparse_vec && ~is_complex
    bag_int_increase = 8;
    bag_ext_increase = 4;
    N_max = Inf;
% elseif isstruct(w)
%     bag_int_increase = 5;
%     bag_ext_increase = 2;
%     N_max = 10;
else
    bag_int_increase = 3;
    bag_ext_increase = 1;
    N_max = 10;
end


N_bag = bag_int_increase;

%atom regeneration
if ~isfield(opt, 'regenerate') 
    opt.regenerate = 1;
end

if ~isfield(opt, 'regen_depth') 
    opt.regen_depth = 0;
end

%visualization
if ~isfield(opt, 'visualize') 
    opt.visualize = 0;
end

if ~isfield(opt, 'visualize_end')
    opt.visualize_end = 0;
end

if ~isfield(opt, 'visualize_delay') 
    opt.visualize_delay = 0;
end

if isfield(opt, 'DG_tol') 
    BM.bagger.DG_tol = opt.DG_tol;
else
    BM.bagger.DG_tol = 1e-4;
end

if ~isfield(opt, 'visualize_delay') 
    opt.norm_tol = 1e-8;
end


will_visualize = opt.visualize || opt.visualize_end || opt.visualize_delay;
first_boundary = 0;


    run_log.error_list  = [];
    run_log.atomic_norm = [];
    run_log.atomic_norm_true = [];
    run_log.cardinality = [];
    run_log.num_attempted = [];
    run_log.num_survived = [];
    run_log.duality_gap = [];
    run_log.time = [];
if  will_visualize
    figure(5);
    clf
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
    % k > 1 && (rank(full([S_new S_bag])) ~= (size(S_new, 2)  + 1))
    if isempty(S_bag)
        
        %check for atomic norm violation
        c = BM.get_c(y);
        anorm_viol = abs(opt.tau*sum(c) - n)>=1e-8;
        
        if opt.regenerate && anorm_viol && BM.update_ext.exterior && ~opt.regen_depth
            %regenerate a new set of atoms matching the current point
            %use a basis pursuit algorithm            
            if isstruct(opt.w)
                
                [BM] = BM.ext_to_full(y);
                y = [1-sum(y); y];
                
                %identify which atoms came from which groups
                %
                %This only works for non-overlapping groups, will need to
                %fix this. Future work?
                Ngroups = length(opt.w.groups);
                group_viol = struct;
                group_viol.group = {};
                group_viol.group_ind = [];
                group_viol.atom_ind= {};
                group_viol.S_rec = {};
                group_viol.c_rec = {};
                group_viol.S_rec_embed = [];
                for g = 1:Ngroups
                    ind_curr = opt.w.groups{g};
                    S_curr = S_new(ind_curr, :);
                    atom_group_curr = find(any(S_curr, 1));
                    atom_group_orig(atom_group_curr) = g;
                    
                    c_curr = c(atom_group_curr);

                    n_ref = opt.tau  * norm(c_curr, 1);                                        
                                        
                    if abs(n_ref - nlist(g)) >= 1e-8
                        %Does this group have an inefficient
                        %representation?
                        weight_curr = opt.w.weights(g);
                        
                        group_viol.group{end+1} = ind_curr;
                        group_viol.group_ind(end+1) =  g;
                        group_viol.atom_ind{end+1} =  atom_group_curr;
                        
                        S_active = S_curr(:, atom_group_curr);
                        x_active = S_active * c_curr;
                        opt_rec = opt;
                        %opt_rec.w = opt.w.weights(g);      
                        opt_rec.w = 1;      
                        opt_rec.num_var = length(x_active);
                        opt_rec.tau = anorm_1d(x_active,  opt_rec.norm_type, opt_rec.w);
                        opt_rec.visualize = 0;
                        opt_rec.visualize_end=0;

                        
                        DG_tol = 1e-8;
                        x_rec = [];
                        
                        while isempty(x_rec) || norm(x_rec - x_active) > 1e-6
                            [x_rec, S_rec, c_rec] = BB_regenerate(x_active, opt_rec, DG_tol);
                            DG_tol = DG_tol * 1e-2;
                        end
                            
                        
                        c_rec = c_rec*weight_curr*opt_rec.tau;
                        S_rec = S_rec/weight_curr;                        
                        
                        S_rec(:, abs(c_rec)<=1e-8) = [];
                        c_rec(abs(c_rec)<=1e-8) = [];
                        
                        S_rec_embed = sparse(opt.num_var, size(S_rec, 2));
                        
                        %embed in larger matrix
                        S_rec_embed(ind_curr, :) = S_rec;
                        group_viol.S_rec{end+1} =  S_rec;
                        group_viol.c_rec{end+1} =  c_rec;
                        group_viol.S_rec_embed = [group_viol.S_rec_embed  S_rec_embed];
                    end
                    
                end
                
                %Now delete the old atoms and add the new atoms back into
                %the system.
                ind_delete = cat(2, group_viol.atom_ind{:});
                c_add = cat(1, group_viol.c_rec{:});
                
                S_rec = group_viol.S_rec_embed;
                BM = BM.delete_indices(ind_delete);
                c_rem = c;
                c_rem(ind_delete) = [];
                y = [c_rem; c_add];
                N_add = length(c_add);
            
            else
                [x_rec, S_rec, c_rec] = BB_regenerate(x, opt);                                                       

                S_rec(:, abs(c_rec)<=1e-8) = [];
                c_rec(abs(c_rec)<=1e-8) = [];
                
                y = c_rec;
                N_add = length(y);
                
                
                %now flush the bash manager
                BM = bash_manager(b_bash, tau, delta, norm_type, w, FCFW);
            end
                                
        
            
                if A_matrix
                    AS_rec = A(S_rec);
                else
                    AS_rec = zeros(length(Ax), N_add);
                    for i = 1:N_add
                        AS_rec(:, i) = A(S_rec(:, i));
                    end
                end
                
                BM = BM.add_atoms(S_rec, AS_rec);
                
                if abs(sum(y)-1) <=  1e-12
                    BM = BM.full_to_ext(c_rec);
                    y = y(2:end);
                end                
            
            [BM, y_new] = BM.bash(y);
            c_new = BM.get_c(y_new);
            
            x_new = BM.get_x(y_new);
            [n_new, nlist_new] = anorm_1d(x_new, opt.norm_type, w);
            Ax_new = BM.get_Ax(y_new);
            S_new = BM.get_S();
            
            error_new = 0.5*norm(Ax_new-b)^2 + 0.5*delta*norm(x_new)^2;
            error_gap = error_old - error_new;
        
            %grad_new = BM.grad(A, y_new);
        
            grad_new = At(Ax_new-b) + delta*x_new;                    
            
        else
            terminate = 1;
            y_new = y;
            error_new = error_old;
            N_survived = 0;
            n_new = n;
            nlist_new = nlist;
        end
        %DG = 0;        
        
        %input is all zeros, or the optimal point is 0
        if k == 1 
            grad_new = grad;
        end
    else
        %DG = real(-grad' * (tau*S_bag(:, 1) - x));
        AS_bag = A(S_bag);
        
        if is_complex
            S_bag = complex_unfold(S_bag, 1);
            AS_bag = complex_unfold(AS_bag, 1);
        end
        
        %BASH
        %Find an optimal loading over current and new atoms
        [BM, y_new] = BM.bash(y, S_bag, AS_bag);
        c_new = BM.get_c(y_new);
        
        on_boundary_new = BM.update_ext.exterior;
        
        %Find statistics of output, including properties of x_new
        
        if first_boundary == 0 && on_boundary_new
            first_boundary = k;
        end
        
        %N_survived = length(c_new) - BM.orig_index;
        %N_dropped  = length(c) - BM.orig_index;
        %fix bag size selection later
        N_survived = size(S_bag, 2);
        N_dropped = 0;
        %update the bag size
        N_bag_next = bag_size_update(N_survived, N_dropped, bag_int_increase, bag_ext_increase, ...
            on_boundary, on_boundary_new, N_max);
        
        x_new = BM.get_x(y_new);
        [n_new, nlist_new] = anorm_1d(x_new, opt.norm_type, w);
        Ax_new = BM.get_Ax(y_new);
        S_new = BM.get_S();
                
        
        if is_complex
            x_new = complex_fold(x_new, 1);
            Ax_new = complex_fold(Ax_new, 1);
            S_new = complex_fold(S_new, 1);
        end
        
        error_new = 0.5*norm(Ax_new-b)^2 + 0.5*delta*norm(x_new)^2;
        error_gap = error_old - error_new;
        
        %grad_new = BM.grad(A, y_new);
        
        grad_new = At(Ax_new-b) + delta*x_new;
        
    end
    
    
    %% logging of the run
    c_new = BM.get_c(y_new);
    run_log.error_list(k)  = error_new;
    run_log.atomic_norm(k) = tau*norm(c_new, 1);
    run_log.cardinality(k) = nnz(c_new);
    run_log.num_attempted(k) = size(S_bag, 2) + length(c);
    run_log.num_survived(k) = N_survived;
    run_log.duality_gap(k) = DG;   
    run_log.time(k) = toc;
    run_log.atomic_norm_true(k) = n_new;
    
    %% visualizations go here
    
    if opt.visualize || (terminate && opt.visualize_end) || (opt.visualize_delay == k)
        %the next project: putting this into a BB_opt.visualizer class
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
            stem(find(~x_new), real(grad_new(x_new == 0)), '.', 'Color', [21 200 225]/255.0)
            stem(find(x_new),  real(grad_new(x_new ~= 0)), '.k')

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
            if is_complex
                x_active = find(x_new);
                for i =  1:length(x_active)
                    ind_curr = x_active(i);
                    x_curr = x_new(ind_curr);
                    plot3([ind_curr, ind_curr], [0, real(x_curr)], [0, imag(x_curr)], 'color', [0.8500    0.3250    0.0980])    
                end
                scatter3(x_active, real(x_new(x_active)), imag(x_new(x_active)), 100, [0.8500    0.3250    0.0980],  '.');
                
                xlabel('$j$', 'Interpreter', 'latex')        
                ylabel('Re $x_j$', 'Interpreter', 'latex')
                zlabel('Im $x_j$', 'Interpreter', 'latex')
                
                pbaspect(aspect_complex)
                view(view_complex)
            else
                stem(x, '.')
                stem(x_new, '.')
                xlabel('$j$', 'Interpreter', 'latex')        
                ylabel('$x_j$', 'Interpreter', 'latex')
            end
        
        hold off
        

        title('Regressor')
        
    %Atomic Norm
        subplot(4, 2, 4)
        hold on
        plot(0:length(run_log.atomic_norm), [atomic_norm_start run_log.atomic_norm])
        plot(0:length(run_log.atomic_norm_true), [atomic_norm_start run_log.atomic_norm_true])
        %stem(x_new, '.')
        if first_boundary
                plot([first_boundary, first_boundary], ylim, ':b')
        end
        plot(xlim, [tau tau], ':k')
        hold off
        
        
        xlabel('iteration')        
        ylabel('$||x||_\mathcal{A}$', 'Interpreter', 'latex')
        title('Atomic Norm')
    
    %Basis Functions
        subplot(4, 2, 5)
        hold on        
        %color_list
        if is_complex
            N_atoms = size(S_new, 2);
            for ka = 1:N_atoms
                c_curr = c_new(ka);
                i_curr = find(S_new(:, ka));
                color_ind = mod(ka-1, length(color_list))+1;
                if norm_type == 1
                    color_curr = 'k';
                else
                    color_curr = color_list(color_ind, :);
                end
                for i = 1:length(i_curr)
                    %too many indices
                    ic = i_curr(i);
                    S_curr = S_new(ic, ka);
                    plot3([ic, ic], tau*c_curr*[0, real(S_curr)], tau*c_curr*[0, imag(S_curr)], 'Color', color_curr)
                end
                scatter3(i_curr, tau*c_curr*real(S_new(i_curr, ka)), tau*c_curr*imag(S_new(i_curr, ka)), 100, color_curr, '.')
            end
            xlabel('$j$', 'Interpreter', 'latex')        
            ylabel('Re $S c_i$', 'Interpreter', 'latex')
            zlabel('Im $S c_i$', 'Interpreter', 'latex')
            
            pbaspect(aspect_complex)
            view(view_complex)
        else
            if BM.bagger.sparse_vec
                basis_func = x_new;
            else
                basis_func = S_new*diag(c_new);
            end
            for i = 1:size(basis_func, 2)
                basis_curr = basis_func(:, i);
                active_curr = find(basis_curr);
                stem(active_curr , basis_curr(active_curr), '.')
            end
                    
            xlabel('$j$', 'Interpreter', 'latex')        
            ylabel('$x_i = S c_i$', 'Interpreter', 'latex')
        end
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
        plot(0:length(run_log.cardinality), [cardinality_start run_log.cardinality])
        plot(run_log.num_attempted, 'xm')
        %stem(x_new, '.')
        if first_boundary
                plot([first_boundary, first_boundary], ylim, ':b')
        end
        hold off
        
        xlabel('iteration')        
        ylabel('$||c||_0$', 'Interpreter', 'latex')
        title('Cardinality')
        
        
    %look at the output
        if opt.visualize || (opt.visualize_delay == k)
             keyboard;
        end
    end
    
    %% prepare for new iteration
    if ~terminate
        k = k+1;
        x = x_new;
        grad = grad_new;
        
        n = n_new;
        nlist = nlist_new;
        
        c = c_new;
        y = y_new;
        error_old = error_new;
        on_boundary = on_boundary_new;
        S_bag_old = S_bag;
        N_bag = N_bag_next;
    end
end

%final output
%package and ship out
x_final = BM.get_x(y);
S_final = BM.get_S();
c_final = BM.get_c(y);

if isfield(opt, 'export_warm_start')
    run_log.warm_start = struct;
    run_log.warm_start.bash_manager = BM;
    run_log.warm_start.y = y;
end

% run_log.warm_start.S = S_final;
% run_log.warm_start.AS =  BM.get_AS();
% run_log.warm_start.K = BM.K;
% run_log.warm_start.rhs = rhs;
% run_log.warm_start.y = y;
% run_log.warm_start.update_ext = BM.update_ext;
% run_log.warm_start.bagger = BM.bagger;

if is_complex
    x_final = complex_fold(x_final, 1);
    S_final = complex_fold(S_final, 1); 
end