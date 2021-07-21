%% Example: MIMO identification
% Let E be prediction error matrix, with Ns rows and ny columns, where Ns =
% number of data samples and ny = number of outputs. Then the minimization
% objective is trace(E*E'*W) where W is a weighting matrix. By default W =
% eye(ny).

function [out, out_random] = exTrace_BB(sys,Ns,SNR,InputBWFraction,opt)
close all
[z,zn,noise,ny,nu,nx] = utGenData(sys,SNR,InputBWFraction,2*Ns);
%{
[ny,nu] = size(sys);
nx = order(sys);
s = sqrt(db2mag(SNR));
u = randn(Ns,nu);
y = lsim(sys,u); % note: zero IC
n = randn(Ns,ny);
yn = y;
for ky = 1:ny
   yn(:,ky) = y(:,ky) + norm(y(:,ky),2)/s*n(:,ky);
end
z = iddata(y,u,1);
zn = iddata(yn,u,1);
%}
ze = zn(1:Ns);
yn = ze.y;
u = ze.u;
y =  z(1:Ns).y;

G = etfe(iddata(y,u));
opt.FreqSample = G.Frequency;
%opt.FreqResponse = permute(G.ResponseData, [3, 2, 1]);
opt.FreqResponse = G.ResponseData;
% 
% TargetCost = norm(yn-y)^2/2;
% fprintf('--------------------------------------------\n')
% fprintf('Actual residue norm: %g\n',TargetCost)
% fprintf('--------------------------------------------\n')

[out, out_random] = atomic_LTI_iteration(y, u opt);

syse = out.sys_out;

e = In.ym-out.y;
fprintf('Generating benchmark results ...')
m=tf(tfest(fft(zn),nx,nx,'feed',1,'ts',zn.Ts,tfestOptions('enforce',1)));
mt=tf(tfest(zn,nx,nx,'feed',1,'ts',zn.Ts,tfestOptions('enforce',1)));
n=ss(ssest(fft(zn),nx,'feed',1,'ts',zn.Ts,ssestOptions('enforce',1)));
nt=ss(ssest(zn,nx,'feed',1,'ts',zn.Ts,ssestOptions('enforce',1)));
nss=ss(ssregest(zn,nx,'feed',1,'ts',zn.Ts));
%fprintf('done.\n')
%fprintf('Reducing model order...')
%sysr=balred(ss(syse),nx,balredOptions('State','truncate'));
z2 = iddata(out.y,u,1);
sysr2=n4sid(z2,nx, 'Feed', 1);
fprintf('done.\n')

out.Results = struct('In',In,'z',z,'zn',zn,...
   'syse',syse,'sysr2', sysr2,...
   'sys_tfest_fd',m,'sys_tfest_td',mt,...
   'sys_ssest_fd',n,'sys_ssest_td',nt,...
   'sys_ssregest_td',nss);

disp('-----------------------------------------')
[~,fit] = compare(z,m,mt,n,nt,syse,sysr2,'init','z');
fprintf('Error norm: %g, Fits (last 2 rows are ours): \n',(norm(e)^2)/2)
disp(cell2mat(fit')')
disp('-----------------------------------------')
end