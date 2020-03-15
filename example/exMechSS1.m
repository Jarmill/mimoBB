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
clear,close all,clc
warning off 'MATLAB:nearlySingularMatrix'
warning off 'MATLAB:singularMatrix'

%{
% enable this when the code is ready
load ASML1Data
Factor = 1e4;
y = StepResp.y*Factor;
G=fselect(G,[1:100:1965,1966:3655]);
G.ResponseData = Factor*G.ResponseData;
G.Frequency = G.Frequency/1000;
z=iddata(y,ones(size(y)),1,'Tstart',0);
zd = z(1:2e4);
%}

rng(2)
sys = rss(4);

%pole(sys)
Ts = 0.04;
t = (0:Ts:25)';
y=step(sys,t);
z = iddata([zeros(10,1);y],[zeros(10,1);ones(size(y))],Ts);
w = logspace(-1,1,50);
G=idfrd(frd(sys,w));
G.Frequency = G.Frequency/25;
z.Ts=1;

%zd=resample(z,1,2);
%zd.Ts = 1;
%opt=ssestOptions;opt.WeightingFilter='inv';opt.Enforce=true;
ny = 1;
nu = 1;
SOLVE = 1;
DRAW = 1;
Nf = numel(G.Frequency);
Ns = size(z,1);
%W = (ones(ny,nu,Nf)./max(1e-6,sqrt(abs(G.ResponseData))))*Ns/Nf;
W = []; %ones(ny,nu,Nf)*Ns/Nf;

opt = sisoAtomOptions;
opt.r1 = 0.9;
opt.phi2 = 0.25*pi;

opt.FreqWeight = W;
%opt.FreqWeight = ones(Ns, ny);
%opt.Compare = 1;
opt.IncludeConstant = 0;
opt.Compare = 0;
opt.tau = 2;
opt.RandomRounds = 20;
opt.ReweightRounds = 20;
opt.NumAtoms = 1000;
opt.NormType = Inf;
opt.FormSystem = true;
if SOLVE
   [out, out_random] = localSolve(z,G,opt); % target cost: 83202.8
   
   if opt.Compare
      utGenAnalysisPlots(out,sys) % quality analysis
   end
end

s=out.sys_modes;
s=cellfun(@(x)ss(x),s,'uni',0);
syse=s{1};
for ct = 2:numel(s)
   syse=syse+s{ct};
end


if DRAW
   figure
   FS = 10;
   clf
   % cm_viridis=viridis(m);
   % colormap(cm_viridis)
   %iopzmap(out.sys_out)
   subplot(1,2,1)
   hold on
   th = linspace(0, 2*pi,  201);
   scatter(real(out_random.poles_active), imag(out_random.poles_active), 200,'x')
   plot(cos(th), sin(th), 'k')
   text(-0.4, 0, sprintf('Order:\nCost:'), 'Fontsize', FS)
   text(0.1, 0, sprintf('%i\n%0.2e', out_random.system_order, out_random.cost), 'Fontsize', FS)
   hold off
   axis square
   box off
   title('Poles before reweighting', 'Fontsize', FS)
   
   xlabel('Re(z)')
   ylabel('Im(z)')
   xticks([-1,-0.5,0,0.5,1])
   yticks([-1,-0.5,0,0.5,1])
   
   
   subplot(1,2,2)
   hold on
   th = linspace(0, 2*pi,  201);
   scatter(real(out.poles_active), imag(out.poles_active), 200, 'x')
   plot(cos(th), sin(th), 'k')
   text(-0.4, 0, sprintf('Order:\nCost:'), 'Fontsize', FS)
   text(0.1, 0, sprintf('%i\n%0.2e', out.system_order, out.cost), 'Fontsize', FS)
   hold off
   axis square
   title('Poles after reweighting', 'Fontsize', FS)
   xlabel('Re(z)')
   ylabel('Im(z)')
   xticks([-1,-0.5,0,0.5,1])
   yticks([-1,-0.5,0,0.5,1])
   
   figure
   bodemag(G,syse)
   
   figure
   compare(z,syse,'init','z')
end

%------------------------------------------------------------------------------------
function [out, out_random] = localSolve(z,G,opt)

nx = 11844; ny = 1; nu = 1;
Ns = size(z,1);
opt.FreqSample = G.Frequency;
opt.FreqResponse = G.ResponseData;
fprintf('Starting Randomization\n')
[ha,p, ~, groups, f, L] = createAtoms(Ns,opt);
W = opt.FreqWeight;
In  = struct('ym',z.y,'u',z.u,'Ts',z.Ts,'ImpRespArray',ha,'PoleArray',p,...
   'PoleGroups', groups, 'FreqRespArray', f, 'FreqWeight', W, ...
   'warm_start', 1, 'ConstWeight', L,...
   'NormType', opt.NormType);

In.ImpRespArray = ha;
In.FreqRespArray = f;
%In.sys_modes = {};
cost_list = [];

if ~opt.RandomRounds && ~opt.ReweightRounds
   opt.FormSystem = 1;
end

out = atomic_LTI_BB(In,opt);
cost_list_random = [cost_list; out.cost];
fprintf('Cost: %0.3e \t Order: %i \t Time: %0.4f \n',out.cost, out.system_order, out.time)

%Multiple rounds of randomizing poles
TRUE_WARM_START = 1;

for i = 1:(opt.RandomRounds)
   %warm_start = out.run_log.warm_start;
   cost_old = out.cost;
   
   active_ind = out.poles_active_ind;
   
   ha_curr =  In.ImpRespArray(:, active_ind);
   f_curr  =  In.FreqRespArray(:, active_ind);
   
   opt.IncludeConstant = 0;
   [ha_new,p_new, scales_new, groups_new, f_new] = createAtoms(Ns,opt);
   
   In.ImpRespArray = [ha_curr ha_new];
   In.FreqRespArray = [f_curr f_new];
   In.PoleGroups = [out.PoleGroups_new (groups_new + max(out.PoleGroups_new))];
   In.PoleArray = [reshape(out.poles_active, 1, []) p_new];
   
   if TRUE_WARM_START
      %prepare warm start for new set of atoms
      active_group_ind_cell = arrayfun(@(i) out.w.groups{i}, out.group_active,   'UniformOutput', false);
      active_group_ind = cat(2, active_group_ind_cell{:});
      warm_start = out.run_log.warm_start;
      BM = warm_start.bash_manager;
      %S_prev = BM.get_S();
      S_prev = BM.S;
      atom_size_max = BM.atom_size_max;
      S_active = S_prev(active_group_ind, :);
      
      S_new = [S_active; sparse(nu*ny*length(p_new), atom_size_max)];
      BM.S = S_new;
      
      
      %this is replicating the anchor atom
      %bug fixing!
      if BM.update_ext.exterior
         S_anc_prev = BM.update_ext.S_anc(active_group_ind);
         S_anc_new = [S_anc_prev; sparse(nu*ny*length(p_new), 1)];
         BM.update_ext.S_anc = S_anc_new;
      end
      
      
      warm_start.bash_manager = BM;
      In.warm_start = warm_start;
   end
   %
   
   if ~opt.ReweightRounds && (i == opt.RandomRounds)
      opt.FormSystem = 1;
   end
   %
   out = atomic_LTI_BB(In,opt);
   %In.sys_modes = out.sys_modes;
   cost_list_random = [cost_list_random; out.cost];
   fprintf('Cost: %0.3e \t (dCost = %0.3e) \t Order: %i \t  Time: %0.4f \n', out.cost, out.cost - cost_old, out.system_order, out.time)
   
end

%In = rmfield(In, 'warm_start');
cost_old = out.cost;
fprintf('Starting Reweighting\n')
out.cost_list_random = cost_list_random;

out_random = out;

cost_list_reweight = [];
%In.warm_start = 1;
%Sparsify the resultant poles through reweighted heuristic
for i = 1:(opt.ReweightRounds)
   
   cost_old = out.cost;
   active_ind = out.poles_active_ind;
   opt.tau = out.tau;
   
   
   if TRUE_WARM_START
      %prepare warm start for reweighted heuristic
      %change to a skewed atomic ball while preserving the same point
      Ngroups = length(out.group_active);
      active_group_ind_cell = arrayfun(@(i) out.w.groups{i}, out.group_active,   'UniformOutput', false);
      active_group_ind = cat(2, active_group_ind_cell{:});
      
      warm_start = out.run_log.warm_start;
      BM = warm_start.bash_manager;
      y = warm_start.y;
      
      was_boundary = BM.update_ext.exterior;
      if BM.update_ext.exterior
         y = [1-sum(y); y];
         BM = BM.ext_to_full();
      end
      
      x_prev = BM.get_x(y);
      S_prev = BM.get_S();
      %S_prev = BM
      %find out which atom came from which group
      atom_group_orig = zeros(Ngroups,1);
      for k = 1:Ngroups
         ind_curr = active_group_ind_cell{k};
         S_curr = S_prev(ind_curr, :);
         atom_group_curr = find(any(S_curr, 1));
         atom_group_orig(atom_group_curr) = k;
         %y_sum(k) = sum(y(atom_group_curr));
      end
      atom_group_num = accumarray(atom_group_orig, 1);
      y_sum = accumarray(atom_group_orig, y );
      
      %now reweight the atoms
      atom_size_max = BM.atom_size_max;
      
      W_ratio = out.PoleGroupWeights_old ./ out.PoleGroupWeights_new;
      W_active  = W_ratio(atom_group_orig)';
      
      
      Natoms = length(y);
      W_diag = sparse(1:Natoms, 1:Natoms, W_active, Natoms, Natoms);
      
      %I think this is correct
      y_w = y ./ W_active;
      
      %Reweight the atoms
      S_prev = BM.S(:, 1:Natoms);
      S_prev_w = S_prev * W_diag;
      %S_new = S_prev_w(active_group_ind, :);
      BM.S(:, 1:Natoms) = S_prev_w;
      BM.S = BM.S(active_group_ind, :);
      
      
      %Reweight the data contributions
      AS_prev = BM.AS(:, 1:Natoms);
      AS_prev_w = AS_prev * W_diag;
      BM.AS(:, 1:Natoms) = AS_prev_w;
      
      %Reweight the kernel and answer vector
      K_w = W_diag * BM.K * W_diag;
      BM.K = K_w;
      
      rhs_w = W_active .* BM.rhs;
      BM.rhs = rhs_w;
      
      if was_boundary
         y_w = y_w(2:end);
         BM  = BM.full_to_ext();
      end
      
      x_w = BM.get_x(y_w);
      grad_w = BM.gradient_full(y_w);
      e_w = BM.get_error(y_w);
      
      %Move onto a stable-optimal loading on the new skewed ball
      %do the bashing inside of optimization
      %BM0 = BM;
      %[BM, y_bash] = BM.bash(y_w);
      
      %x_bash = BM.get_x(y_bash);
      %grad_e = BM.gradient_full(y_bash) ;
      %e_bash = BM.get_error(y_bash);
      %this is replicating the anchor atom
      
      warm_start.bash_manager = BM;
      warm_start.y = y_w;
      In.warm_start = warm_start;
   end
   
   ha_curr =  In.ImpRespArray(:, active_ind);
   f_curr  =  In.FreqRespArray(:, active_ind);
   
   In.ImpRespArray = ha_curr;
   In.FreqRespArray = f_curr;
   
   In.h0 = out.h;
   In.PoleGroupWeights = out.PoleGroupWeights_new;
   In.PoleGroups = out.PoleGroups_new;
   In.PoleArray = out.poles_active;
   
   %In.PoleGroupWeights = out.PoleGroupWeights_new_all;
   
   
   if (i == opt.ReweightRounds)
      opt.FormSystem = 1;
   end
   
   out = atomic_LTI_BB(In,opt);
   %In.sys_modes = out.sys_modes;
   fprintf('Cost: %0.3e \t (dCost = %0.3e) \t Order: %i \t  Time: %0.4f \n', out.cost, out.cost - cost_old, out.system_order, out.time)
   
   cost_list_reweight = [cost_list_reweight; out.cost];
   
   if abs(out.cost - cost_old) <= opt.ReweightTol
      break
   end
   
   cost_old = out.cost;
end

out.cost_list_reweight = cost_list_reweight;

%Extract system

c=out.Coeff{1};
p=out.poles_active;


if opt.Compare
   % Num = out.h';
   % Num = cellfun(@(x)x.',Num,'uni',0);
   % syse = tf(Num,num2cell(ones(ny,nu)),1,'var','z^-1');
   
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
end