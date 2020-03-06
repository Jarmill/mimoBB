function [a, n] = LMO_1d(G, norm_type, w)
%LMO_1D: Evaluates the linear minimization oracle to find the best atom to
%add, max <G, a> for a in the set A. Returns a single atom.
%
%Input:
%   G:          Direction to try and match (full or projected -gradient)
%   norm_type:  What norm to use
%   w:          Additional information in the norm
%
%Output:
%   a:          Best atom chosen by LMO.
%   n:          Output norm

%Atoms a and the 'negative gradient' G come from the same Banach space, which means
%they have equivalent form (for vectors at least)
G0 = G;

if isnumeric(norm_type)
    %numerical norm_type means that an Lp norm is used, where p>=1
    
    %For the Lp norms, w is stored with reweighted heuristic values
    %(diagonal scaling). Divide by these weights at the end.
    
    %sum of norms regularization
    %w is a struct with fields 'weights' and 'groups'
    if isstruct(w)
        %sum of norms
        num_groups = length(w.groups);
        norm_list = zeros(num_groups, 1);
        
        %poll all groups for their weighted norm
        n = 0;
        for i = 1:num_groups
            g_curr = w.groups{i};
            w_curr = w.weights(i);
            norm_curr = norm(G(g_curr), norm_type);
            norm_list(i) = norm_curr/w_curr;
            n = n + norm_curr*w_curr;
        end
        
        %find the group with the highest weighted norm
        %can be done without storage, make more efficient later. single
        %O(n) call.
        [~, i_max] = max(norm_list);
        
        %find the best atom in the best group
        g_max = w.groups{i_max};
        w_max = w.weights(i_max);
        
        G_max = G(g_max);
        a_max = LMO_1d(G_max, norm_type, 1);
        
        a = sparse(size(G, 1), 1);
        a(g_max) = a_max / w_max;
    else
        %standard regularization
        weighted = ~isequal(w, 1);    
        n = norm(w.*G, norm_type);
        G = G./w;
        if norm_type == Inf
            %L infinity (corners of ball)                        
            if isreal(G)
                a = sign(G);                
            else
                %normalized phase
                a = G./abs(G);
            end
            
            %think this causes degeneracy?
            %a(abs(G)<=1e-6) = 0;
        elseif norm_type == 1
            %L1 norm (sparse vectors)
            a = sparse(size(G, 1), 1);
            [~, i] = max(abs(G));
            if isreal(G)
                a(i) = sign(G(i));
            else
                %normalized phase
                a(i) = G(i)/abs(G(i));
            end
            
        else
            %Lp norm, 1<p<Inf
            p = norm_type;

            %Ugly code here comes from the subgradient of the dual norm
            
            ag = abs(G);

            if isreal(G)
                sg = sign(G);
            else
                sg = G./abs(G);
            end
            
            %only the magnitude changes, not the phase
            %this doesn't make sense, but roll with it.
            num_a = ag.^(1/(p-1));
            den_a = norm(num_a, p); %normalizes the atom to unit L_p norm

            a = sg .* num_a / den_a;
            
        end

        %if weighted
        a = a./w;
    end    
else
    if strcmp(norm_type, 'simp')
        %simplex constraint
        %not centrally symmetric, a gauge
        n = sum(w.*G);
        G = G./w;
        ip = find(G>0);
        a = sparse(size(G, 1), 1);
        [~, i] = max(G(ip));       
        a(ip(i)) = sign(G(ip(i)));
        a = a ./ w;
        
    elseif strcmp(norm_type, 'pos')
        %positive orthant intersected with L infinity ball
        %also a gauge, can be described by points of simplex plus a point
        %of all 1's to form the corner
        n = max(w.*G);
        G = G./w;
        ip = find(G>0);
        a = sparse(size(G, 1), 1);
        %[~, i] = max(G(ip));       
        a(G>0) = 1;
        a = a./ w;
        
    elseif strcmp(norm_type, 'chain')
        %k-chain lasso
        %sum of L2 norm penalty, each group is a chain k elements long
        %credit to Marina Vinyes and Guillaume Obozinski, 2016
        
        %Find the w-long sequence of elements (chain) with maximal L2 norm.
        %w is the elements in the chain
        chain = ones(w, 1);
        
        %this convolution trick only works on an L2 norm penalty
        z = conv(G.^2, chain, 'valid');
        
        [n, i] = max(z);
        n = sqrt(n); %maximal norm (dual norm)
        
        window = i:(i+w-1);
        a_chain = G(window);
        a_chain_norm = a_chain ./ norm(a_chain);
        
        a = sparse(window, chain, a_chain_norm, length(G), 1);

    elseif strcmp(norm_type, 'OWL')
        %OWL norm with weights w
        %is centrally symmetric
        %Ordered Weighted L1 norm
        
        %modify the weights w
        ws = cumsum(w);         %weights are a weakly decreasing sequence
        wr = 1./ws;             %reciprocol of the cumulative sum, is the weights 
                                %on each atom for a set of coordinates chosen

        %find the best atom in the set

        %sort indices by highest absolute value of the gradient
        %may be the dominant cost
        [G_abs_sort, i_max] = sort(abs(G), 'descend');

        %get the cumulative sums of the greatest absolute values
        G_abs_sum = cumsum(G_abs_sort);

        %compare the weighted sum of using k components. The more components
        %used, the higher the cumulative sum, but the lower the associated
        %atomic weight (sum of reciprocals of weights w to get to that point)
        
        %Use the linear minimization oracle to decide (subgradient of dual OWL
        %norm)
        LMO_OWL = wr .* G_abs_sum;
        [n, N_components] = max(LMO_OWL); %how many components to add
        ind_added = i_max(1:N_components); %the indices corresponding to added components


        %the final atom is the signs of G's components weighted by the atomic
        %weight of adding that many components. Kind of confusing, but easier
        %than the OWL paper's hadamard product and permutation and orderings.
        %a = sparse(length(x), 1);
        a = zeros(size(G));
        if isreal(G)
            a(ind_added) = sign(G(ind_added)) * wr(N_components);
        else
            a(ind_added) = G(ind_added)./abs(G(ind_added)) * wr(N_components);
        end
        
    elseif strcmp(norm_type, 'poles')    
        %random sample of unit circle in z domain
        %possible sector conditions as well
        %   annulus between radii rho1 and rho2
        %   maximum angle max_angle
        if ~isfield(w, 'n')
            w.n = 100;
        end
        if ~isfield(w, 'rho1')
            w.rho1 = 0;
        end
        if ~isfield(w, 'rho2')
            w.rho2 = 1;
        end
        if ~isfield(w, 'max_angle')
            w.max_angle = pi;
        end
        
        %sample the poles
        p = uniform_over_ring_sector(w.rho1, w.rho2, w.n, w.max_angle);
        [a,p_best,val_best] = best_atom(G,p);
    else
        %point cloud
        %not centrally symmetric, the atomic norm a gauge defined by how
        %much the convex hull of the point set needs to be stretched to
        %contain x
        %
        %w contains the list of points
        
        LMO_point = G'*w;
        [n, i_max] = max(LMO_point);
        a = w(:, i_max);
    end

end
        