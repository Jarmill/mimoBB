% Script to find the best atom given a gradient vector and a set of 
% poles
% Author : Burak Yilmaz
% Last Update : 10/27/2014

% cumprod is used to raise the atom to powers. Using logarithms 
% is faster for N ~ 2000 or larger. 

function [a,p_best,val_best] = best_atom(grad_f,p)

    N = length(grad_f);
%    scales = get_scalesBen(p,N);
    [scale_r,scale_c] = get_scales(p,N); 
   % scale_r = scales; 
   % scale_c = scales;
    grad_f = grad_f(:)';
    k = length(p);   % Number of poles checked
    index_N = ones(1,N-2);
    A = [zeros(1,k);ones(1,k);cumprod(p(index_N,:))]; % Atomic responses are 
    % simply real(A) for set A1 and imag(A) for set A2, before applying scales.
    
    % Form the dot products
    temp2 = grad_f*A;
    
    prod_r = scale_r.*real(temp2);
    prod_c = scale_c.*imag(temp2);
    [val_r,Ir] = max(abs(prod_r));
    [val_c,Ic] = max(abs(prod_c));
    
       
    if val_r>val_c
        p_best = p(Ir);
        a = -sign(prod_r(Ir))*(2*scale_r(Ir))*real(A(:,Ir));  
        val_best = val_r;
    else
        p_best = p(Ic);
        a = -sign(prod_c(Ic))*(2*scale_c(Ic))*imag(A(:,Ic));  
        val_best = val_c;
    end
   
  
end