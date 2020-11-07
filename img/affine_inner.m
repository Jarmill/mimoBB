function [x_aff] = affine_inner(s, tau, data, ref, use_kernel, on_boundary, delta, w)
%AFFINE_INNER
%finds the optimal point (minimizes function, x_aff) on affine subspace 
%(face) that contains the point x.

%   x:      current point. Optimization will occur over all nonzero coordinates of x
%   tau:    Maximum L1-norm. If norm(x, 1) = tau, on_boundary = 1
%   data:      either K = A'A + delta*I (kernel) or A (non-kernel)
%   ref:    either A'b (kernel) or b (non-kernel)
%   use_kernel: whether to use kernelized formulation
%   on_boundary: whether to use an interior cutting plane step (0) or a slant/facial step (1)
%   w: weights for the reweighted heuristic

%if norm(x_aff, 1) = tau (still on same face), then all nonzero coordinates
%that have been chosen have the correct sign.

%if norm(x_aff, 1) > tau (on affine subspace, but outside L1-ball of radius
%tau), then some of the dimensions that have been chosen are invalid,
%meaning their sign is bad or should be 0. 

%solves a linear system r-1 x r-1, where nnz(x) = r. Can be expensive 

%if nargin < 8
%    w = 1;
%end

s_ind = find(s);

%the good stuff. affine steps. Expected to yield massive gains in
%performance.

%first find out B: the matrix of directions that define the affine subspace
%for reweighted heuristic, is actually W^-1 * B, done automatically
if on_boundary
    %on boundary, so a slant step (swapping between coordinates) is
    %necessary.
    %First find the affine subspace spanned by the corners.

    %not sure if I can vectorize this vs. put it in a for loop
    ise = s_ind(end);

    %reference point on the affine subspace. could be x, but is more
    %efficient sparsitywise to choose one of the corners.
    p = sparse(size(s, 1), 1);
    p(ise) = tau*sign(s(ise));
    
    i_active = s_ind(1:end-1);
    s_active = s(i_active);
    w_active = w(i_active);
    
    %directions that define the affine subspace
    Bnum = length(s_ind)-1; %number of corners supporting affine subspace

    %B = [tau s1 e1 - tau sn en; tau s2 e2 - tau sn en;... tau sn-1 en-1 - tau sn en]
    p_values_w = ones(Bnum, 1) * -(tau*s(ise))/w(ise); %from reference value
    b_values_w = tau * s_active ./ w_active; %other non-reference corners
    p_values = ones(Bnum, 1) * -(tau*s(ise)); %from reference value
    b_values = tau * s_active; %other non-reference corners

    i = [ones(1, Bnum)*ise s_ind(1:end-1)']';
    j = [1:Bnum 1:Bnum]';
    v = [p_values; b_values];
    vw = [p_values_w; b_values_w];

else
    %on the interior, so there is full freedom of travel.
    
    %in this code, are optimizing over all previously chosen dimensions, so
    %the reference point is always zero.
    
    %if a limited selection of axes were being optimized over, p would be
    %filled with the value of the constant axes
    %Example: if x = [x1; x2; x3; x4], and optimization would occur over 
    % dimensions 1 and 3, then p = [0; x2; 0; x4]
    p = zeros(size(s));
    
    %full freedom of motion here
    i_active = s_ind;
    s_active = s(i_active);
    w_active = w(i_active);
    
    
    Bnum = length(s_ind);
    i = i_active;
    j = (1:Bnum)';
    %v = sign(x(x_ind));
    %v = ones(size(j))./w(s_ind);
    v = tau * s_active; %actually equivalent to 1   
    vw = v ./ w_active;
end

B = sparse(i, j, v, size(s, 1), Bnum);
Bw = sparse(i, j, vw, size(s, 1), Bnum);

%is an absolute nightmare, need to clean up and document.

%This algorithm is solving the normal equations, which is the current
%bottleneck of processing. Multiplying AB'*AB is expensive, but at least
%it is smaller than A'*A. Maybe a QR type of method will work or help out?

%B is tall and skinny, so AB'*AB is small, while Q and R are large. QR
%algorithms are therefore inefficient. It is entirely possible that the
%normal equations are the way to go for this problem. Since sparse
%solutions are desired, it is reasonable to assume that B will remain
%extremely tall and skinny for all future applications.

if use_kernel
    %process input
    K = data;
    Atb = ref;
    
    %find linear system
    K_aff = B'*(K*B);    
    rhs = B'*(Atb - K*p);
else
    A = data;
    b = ref;
    
    %find linear system
    AB = A*B;
    %K_aff = (AB'*AB) + delta*(B'*B);    
    
    K_aff_left = AB'*AB;
    K_aff_right = delta*(Bw'*Bw);
    K_aff = K_aff_left + K_aff_right;
    %K_aff is a dense, SPD matrix
       
    rhs = AB'*(b-A*p) - delta*Bw'*(p./w);
end

%basic formula:
%x_aff = B ((B' K B) \ B'(A'b - K p)) + p

%the bottleneck is in computing K_aff, not in finding the solution.
best_y = (K_aff \ rhs);
%best_y = pcg(K_aff, rhs, 1e-8, 50);


%this is the point x that minimizes the elastic net function on the
%affine subspace spanned by its residing face. 
x_aff = sparse(B*best_y + p);

%outer code will figure out if x_aff is a valid point on the face


end

