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


fprintf('Starting Randomization\n')
pp = pole(sys);
[ha,p, scales, groups, f, L] = createAtoms(Ns,opt);
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
         figure('pos',pos{ky,ku},'Name',fprintf('IO(%d,%d)',ky,ku))
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
'TargetCost',TargetCost,'pH',{pH}, 'warm_start', 1, 'ConstWeight', L,...
'NormType', opt.NormType);

In.ImpRespArray = ha;
In.FreqRespArray = f;


cost_list = [];

if ~opt.RandomRounds && ~opt.ReweightRounds
    opt.FormSystem = 1;
end

out = atomic_LTI_BB(In,opt); 
cost_list_random = [cost_list; out.cost];
fprintf('Cost: %0.3e \t Time: %0.4f \n',out.cost, out.time)

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
    cost_list_random = [cost_list_random; out.cost];
    fprintf('Cost: %0.3e \t (dCost = %0.3e)  Time: %0.4f \n', out.cost, out.cost - cost_old, out.time)
    
end

%In = rmfield(In, 'warm_start');
cost_old = out.cost;
fprintf('Starting Reweighting\n')
out.cost_list_random = cost_list_random;

cost_list_reweight = [];
%In.warm_start = 1;
%Sparsify the resultant poles through reweighted heuristic
for i = 1:(opt.ReweightRounds)   

    cost_old = out.cost;
    active_ind = out.poles_active_ind;    
    
    
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
    fprintf('Cost: %0.3e \t (dCost = %0.3e)  Time: %0.4f \n', out.cost, out.cost - cost_old, out.time) 
    cost_old = out.cost;
    cost_list_reweight = [cost_list_reweight; out.cost];
end

out.cost_list_reweight = cost_list_reweight;

%Extract system



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