%% Example: MIMO identification
% Let E be prediction error matrix, with Ns rows and ny columns, where Ns =
% number of data samples and ny = number of outputs. Then the minimization
% objective is trace(E*E'*W) where W is a weighting matrix. By default W =
% eye(ny).

function out = exTrace_BB(sys,Ns,SNR,InputBWFraction,opt)
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

TargetCost = norm(yn-y)^2/2;
fprintf('--------------------------------------------\n')
fprintf('Actual residue norm: %g\n',TargetCost)
fprintf('--------------------------------------------\n')

pp = pole(sys);
[ha,p, scales, groups, f] = createAtoms(Ns,opt);
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


%Need to cleverly design weighting functions (to do)

W = opt.FreqWeight;

In = struct('ym',yn,'u',u,'Ts',1,'ImpRespArray',ha,'PoleArray',p,...
   'PoleGroups', groups, 'FreqRespArray', f, 'FreqWeight', W, ...
'TargetCost',TargetCost,'pH',{pH});

p_curr  = p;
f_curr  = f;
ha_curr = ha;
groups_curr = groups;

for i = 1:(opt.ReweightRounds + 1)
    
    In.ImpRespArray = ha_curr;
    In.FreqRespArray = f_curr;
    
    out = atomic_LTI_BB(In,opt); 
    In.h0 = out.h;
    In.PoleGroupWeights = out.PoleGroupWeights_new;
    In.PoleGroups = out.PoleGroups_new;
    In.PoleArray = out.poles_active;
    active_ind = out.poles_active_ind;    
    
    ha_curr =  ha_curr(:, active_ind);
    f_curr  =  f_curr(:, active_ind);        
end


if opt.Compare
Num = out.h';
Num = cellfun(@(x)x.',Num,'uni',0);
syse = tf(Num,num2cell(ones(ny,nu)),1,'var','z^-1');

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