function [n, norm_list] = anorm_1d(G, norm_type, w)
%evaluate atomic norm at a point G
%LMO performs the dual norm
if nargin < 3
    w = 1;
end

if isstruct(w)
    %Ngroups = length(w.groups);
    %norm_list = zeros(Ngroups, 1);
%     for i = 1:Ngroups                
%         g_curr = w.groups{i};
%         w_curr = w.weights(i);
%         norm_curr = anorm_1d(G(g_curr), norm_type);
%         norm_list(i) = norm_curr*w_curr;
%     end
    norm_list_0 = cellfun(@(wg) anorm_1d(G(wg) ,norm_type), w.groups);
    if size(norm_list_0, 1) ~= size(w.weights, 1)
        norm_list = norm_list_0 .* w.weights';
    else
        norm_list = norm_list_0 .* w.weights;
    end    
    n = sum(norm_list);
else
    if isnumeric(norm_type)
        n = norm(G.*w, norm_type);
    elseif strcmp(norm_type, 'simp')
        n = sum(G.*w);
    elseif strcmp(norm_type, 'pos')
        n = max(w.*G);
        n = max(n, 0)
    elseif strcmp(norm_type, 'chain')
        chain = ones(w, 1);
        z = conv(G.^2, chain, 'valid');
        
        [n, i] = max(z);
        n = sqrt(n); %maximal norm (dual norm)
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
    else
        %point cloud
        %w contains the list of points
        LMO_point = G'*w;
        [n, i_max] = max(LMO_point);
    end
    norm_list =  n;
end