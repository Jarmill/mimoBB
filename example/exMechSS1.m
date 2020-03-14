% Practical example: Pascal

%{
[M,C,K,B] = getMECHSS('piston');
rng(0), F = sprand(1,size(K,1),.01); 
sys1 = mechss(M,C,K,B,F);
t = linspace(0,0.03,5000);
step(sys1,t);
%}

%{
[M,C,K,B,F] = getMECHSS('ASML1');
sys2 = mechss(M,C,K,B,F);
w = logspace(1,7,10000);
sigma(sys2,w); 
%}
