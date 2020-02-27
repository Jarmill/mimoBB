%w = 3;

%G = randn(10, 1);

%a = LMO_1d(G, 'chain', w);
%[alist, N_added] = bag_chain(-G, 0*ones(size(G)), 3, 4, w);


rng(0)
n=300;
p=1000;
% n=600;
% p=2000;

k=8;
rho=5;

% k=1;
% rho = 100;

sigma=0.1;

%% synthetic data
% Random design matrix D, paramter vector x, output y
m=10;
x=[2*ones(m,1);zeros(p-m,1)];

U=randn(n,n);
[U,R]=qr(U);
V=randn(p,p);
[V,R]=qr(V);
S=zeros(n,p);
s=0.95.^(0:1:(min(n,p)-1));
S(1:min(n,p),1:min(n,p))=diag(s);
D=sqrt(n)*U*S*V';
 
y=D*x+sigma*randn(n,1);
[x_bb, S_bb, c_bb, run_log] = BB_1d(D, y, rho, 0, 'chain', k);