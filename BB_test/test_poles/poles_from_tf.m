function [c, p] =  poles_from_tf(sys)
%from system, generates coefficients and poles
%these coefficients are real (sin/cos), not complex conjugate residues

%sys: tf (b/a) representation of the system
b = sys.Numerator{1};
a = sys.Denominator{1};
[r, p, ~] = residue(b, a);
c = r;
%deal with complex poles/complex residues
%as per convention, imag(p)>0: cos, imag(p)<0: sin, imag(p)=0: exp
c(imag(c) > 0) = 2*real(c(imag(c) > 0));
c(imag(c) < 0) = 2*imag(c(imag(c) < 0));


end