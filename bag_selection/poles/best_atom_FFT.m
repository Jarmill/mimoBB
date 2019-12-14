% Script to find the best atom given a gradient vector and a set of 
% poles
% Author : Burak Yilmaz
% Last Update : 11/03/2014

% cumprod is used to raise the atom to powers. Using logarithms 
% is faster for N ~ 2000 or larger. 

function [a,p_best,a_pred,val_best] = best_atom_FFT(grad_f,p,Nf,fft_mat,Nv)

    %Nf = 2^3;   % Nf point FFT . 
    N = length(grad_f);
    index_N  = ones(1,N-2);
    index_Nf = ones(1,Nf);
    temp = p(index_Nf,:).*fft_mat; 
    k = length(p); 
    
    scale_r = zeros(Nf,k);
    scale_c = zeros(Nf,k);
   [scale_r(:),scale_c(:)] = get_scales(temp(:),N); 
%     if norm(scale_r-scale_r1)>1e-5
%         norm(scale_r-scale_r1)
%     end
%     if norm(scale_c-scale_c1)>1e-5
%         norm(scale_c-scale_c1)
%     end
%     
    grad_f = grad_f(:);
          
    A = [2*ones(1,k);2*cumprod(p(index_N,:))];   
    F = grad_f(2:end,ones(1,k));
    temp2 = (fft(F.*A,Nf));  

    prod_r = scale_r.*real(temp2);
    prod_c = scale_c.*imag(temp2);
    [val_r,Ir] = max(abs(prod_r));
    [val_c,Ic] = max(abs(prod_c));
    
    [val,I] = max([val_r;val_c],[],2);
    
    if val(1)>val(2)
        p_best = p(I(1))*fft_mat(Ir(I(1)),I(1));
        dummy2 = [1;cumprod(p_best(index_N,:))];   
        a = -sign(prod_r(Ir(I(1)),I(1)))*(2*scale_r(Ir(I(1)),I(1)))*real(dummy2);  
        val_best = val(1);
        dummy3(1:N+Nv-2,1) = p_best;
        comp_vec3 = cumprod(dummy3,1);    
        a_pred = [0;1;real(comp_vec3)]*(-sign(prod_r(Ir(I(1)),I(1)))*(2*scale_r(Ir(I(1)),I(1)))); 
        
    else
        p_best = p(I(2))*fft_mat(Ic(I(2)),I(2));
        dummy2 = [1;cumprod(p_best(index_N,:))];   
        a = -sign(prod_c(Ic(I(2)),I(2)))*(2*scale_c(Ic(I(2)),I(2)))*imag(dummy2);  
        val_best = val(2); 
     
        dummy3(1:N+Nv-2,1) = p_best;
        comp_vec3 = cumprod(dummy3,1);    
        a_pred = [0;0;imag(comp_vec3)]*(-sign(prod_c(Ic(I(2)),I(2)))*(2*scale_c(Ic(I(2)),I(2)))); 
    end
    a = [0;a];
%     test(:,:,1) = scale_r.*real(prod);
%     test(:,:,2) = scale_c.*imag(prod);
%    
%     [val,I] = max(abs(test(:)));
%     %val*2
%     [i,j,k] = ind2sub(size(test),I);
%     p_best = p(j)*fft_mat(i,1);
%     dummy2 = [1;cumprod(p_best(index_N,:))];   
%     if k==1
%         a = -sign(real(prod(i,j)))*(2*scale_r(i,j))*real(dummy2);  
%     else
%         a = -sign(imag(prod(i,j)))*(2*scale_c(i,j))*imag(dummy2);
%     end
 end