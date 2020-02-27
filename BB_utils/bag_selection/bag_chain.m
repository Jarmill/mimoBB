function [atom_bag_list, N_added] = bag_chain(grad, x, N_bag, tau, w)
%BAG_CHAIN k-chain group lasso
%sum of L2 norm penalty, each group is a chain k elements long
%credit to Marina Vinyes and Guillaume Obozinski, 2016

    %Find the w-long sequence of elements (chain) with maximal L2 norm.
    %w is the elements in the chain
    G = -grad;
    chain = ones(w, 1);

    %this convolution trick only works on an L2 norm penalty
    z = conv(G.^2, chain, 'valid');

    [v, i] = maxk(z, N_bag);
    v = sqrt(v); %maximal norm (dual norm)

    %fill up the atom bag list
    atom_i = zeros(N_bag*w, 1);
    atom_j = zeros(N_bag*w, 1);
    atom_v = zeros(N_bag*w, 1);
    for k = 1:length(i)
        i_curr = i(k);
        window = i_curr:(i_curr+w-1);
        a_chain = G(window);
        a_chain_norm = a_chain ./ norm(a_chain);
        
        atom_i(((k-1)*w + 1): k*w) = window;
        atom_j(((k-1)*w + 1): k*w) = k*chain;
        atom_v(((k-1)*w + 1): k*w) = a_chain_norm;
        %atom_i = [atom_i; window];
        %atom_j = [atom_j; k*chain];
        %atom_v = a_chain_norm;
        
        %atom_bag_list_screen(window, k*chain) = a_chain_norm;
    end
    %atom_bag_list_screen = sparse([], [], [], length(G), N_bag);
    atom_bag_list_screen = sparse(atom_i, atom_j, atom_v, length(G), N_bag);
    
    %find which atoms are feasible
    %dir_travel = tau*atom_bag_list_screen - x;
    %DG = real(dir_travel'*(-grad));

    DG = grad'*x + tau*v;
    survive_index = find(DG > 1e-5);
    N_added = length(survive_index);
    atom_bag_list = atom_bag_list_screen(:, survive_index);
end

