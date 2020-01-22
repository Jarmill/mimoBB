%% Example: MIMO identification
% Let E be prediction error matrix, with Ns rows and ny columns, where Ns =
% number of data samples and ny = number of outputs. Then the minimization
% objective is trace(E*E'*W) where W is a weighting matrix. By default W =
% eye(ny).

function out = exTrace_BB(sys,Ns,SNR,opt)
close all
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

TargetCost = norm(yn-y)^2/2;
fprintf('--------------------------------------------\n')
fprintf('Actual residue norm: %g\n',TargetCost)
fprintf('--------------------------------------------\n')

pp = pole(sys);
[ha,p, scales, groups] = createAtoms(Ns,opt);
pH = cell(ny,nu);
if opt.ShowProgressPlot   
   pos = cell(ny,nu);
   pos{1,1} = [414,667,377,304];
   pos{1,2} = [798,666,363,304];
   pos{1,3} = [1168,666,409,298];
   pos{2,1} = [410,281,393,295];
   pos{2,2} = [807,281,360,296];
   pos{2,3} = [1173,283,406,288];
   for ky = 1:ny
      for ku = 1:nu
         figure('pos',pos{ky,ku},'Name',sprintf('IO(%d,%d)',ky,ku))
         hh = plot(real(pp),imag(pp),'rx',NaN,NaN,'go');
         pH{ky,ku} = hh(2);
         shg
      end
   end
end
In = struct('ym',yn,'u',u,'Ts',1,'ImpRespArray',ha,'PoleArray',p,...
   'PoleGroups', groups, 'TargetCost',TargetCost,'pH',{pH});

out = atomic_LTI_BB(In,opt); In.h0 = out.h;

%I have absolutely no idea what these below operations are doing.
%Task for another day.
%Num = permute(mat2cell(out.h,Ns,ones(1,ny),ones(1,nu)),[2 3 1]);
Num = out.h';
Num = cellfun(@(x)x.',Num,'uni',0);
syse = tf(Num,num2cell(ones(ny,nu)),1,'var','z^-1');

e = In.ym-out.y;
fprintf('Generating benchmark results ...')
m=tf(tfest(fft(zn),nx,nx,'feed',1,'ts',zn.Ts,tfestOptions('enforce',1)));
mt=tf(tfest(zn,nx,nx,'feed',1,'ts',zn.Ts,tfestOptions('enforce',1)));
n=ss(ssest(fft(zn),nx,'feed',1,'ts',zn.Ts,ssestOptions('enforce',1)));
nt=ss(ssest(zn,nx,'feed',1,'ts',zn.Ts,ssestOptions('enforce',1)));
fprintf('done.\n')
fprintf('Reducing model order...')
%sysr=balred(ss(syse),nx,balredOptions('State','truncate'));

z2 = iddata(out.y,z.u,1);
sysr2=n4sid(z2,nx, 'Feed', 1);
fprintf('done.\n')

out.Results = struct('In',In,'z',z,'zn',zn,...
   'syse',syse,'sysr2', sysr2,...
   'sys_tfest_fd',m,'sys_tfest_td',mt,...
   'sys_ssest_fd',n,'sys_ssest_td',nt);

disp('-----------------------------------------')
[~,fit] = compare(z,m,mt,n,nt,syse,sysr2,'init','z');
fprintf('Error norm: %g, Fits (last 2 rows are ours): \n',(norm(e)^2)/2)
disp(cell2mat(fit')')
disp('-----------------------------------------')
