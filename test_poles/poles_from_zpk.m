function [c, p] =  poles_from_zpk(sys)
%from system, generates coefficients and poles
%needs to find residues of expansion first 

%sys: zpk representation of the system
sys_tf = tf(sys);

[c, p] = poles_from_tf(sys_tf);

end