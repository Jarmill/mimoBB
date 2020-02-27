function [B, p] = affine_ext_representation(s, tau, on_boundary)
%B, p representation of affine subspace on face of simplex/L1-ball.
%c = By+p,  y are barycentric (homogoneous) coordinates of points in c.
%Input:
%   s:      signs of coefficients under consideration
%   tau:    radius of L1/ball or height of simplex
%   on_boundary: whether this representation is on interior (free space) or
%           exterior (constrained to a face).

%s = sign(c);

s_ind = find(s);

if on_boundary
    %on boundary, so a facial/exterior step (swapping between coordinates) is
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
    
    %directions that define the affine subspace
    Bnum = length(s_ind)-1; %number of corners supporting affine subspace

    %B = [tau s1 e1 - tau sn en; tau s2 e2 - tau sn en;... tau sn-1 en-1 - tau sn en]
    p_values = ones(Bnum, 1) * -(tau*s(ise)); %from reference value
    b_values = tau * s_active; %other non-reference corners

    i = [ones(1, Bnum)*ise s_ind(1:end-1)']';
    j = [1:Bnum 1:Bnum]';
    v = [p_values; b_values];

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
    
    Bnum = length(s_ind);
    i = i_active;
    j = (1:Bnum);
    %v = sign(x(x_ind));
    %v = ones(size(j))./w(s_ind);
    v = tau * s_active; %actually equivalent to 1   
end

%parameters that define affine subspace
B = sparse(i, j, v, length(s), Bnum);
%p = sparse(size(s, 1), 1);
%p(ise) = sign(s(ise));