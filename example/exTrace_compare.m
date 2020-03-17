function [out, out_fw] = exTrace_compare(sys,Ns,SNR,InputBWFraction,opt)
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

[ha,p, scales, groups, f, L] = createAtoms(Ns,opt);


%Need to cleverly design weighting functions (to do)

W = opt.FreqWeight;

In  = struct('ym',yn,'u',u,'Ts',1,'ImpRespArray',ha,'PoleArray',p,...
     'PoleGroups', groups, 'FreqRespArray', f, 'FreqWeight', W, ...
     'warm_start', 1, 'ConstWeight', L,...
     'NormType', opt.NormType);

In.ImpRespArray = ha;
In.FreqRespArray = f;

if ~opt.RandomRounds && ~opt.ReweightRounds
    opt.FormSystem = 1;
end

out = atomic_LTI_BB(In,opt); 

%opt_fw = opt;
%opt_fw.FCFW = 0;
%out_fw = atomic_LTI_BB(In,opt_fw);

opt_fw = atomic_LTI(In, opt);

end