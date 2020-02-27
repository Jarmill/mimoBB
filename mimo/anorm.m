function a = anorm(h)
% compute nuclear norm of a Hankel matrix composed from an impulse response
% vector h. If h is a matrix, the norm is computed column-wise. 

a = 1./getSSScales(h);
