%% MIMO IDENTIFICATION: L2, TIME and FREQUENCY DOMAINS
% Minimize prediction error trace L2 norm for MIMO LTI identification.
% We consider a partial fraction structure G(z) = \sum_i R_i/(z-p_i), where
% R is residue matrix of size [ny, nu]. For mechanical structures, R is a
% dense matrix. But in general, we impose low cardinality of vec(R).
%
%Fully Corrective Bag-and-bash implementation
%
%Use sum-of-norms regularization for MIMO system


% Input fields:
% u     = Input Toeplitz matrix
% ym     = Measured signal (Output)
% tau    = Bound on atomic norm of the solution
% t_max  = Maximum iteration number
% p      = Random poles
% m_pole = Number of poles to try @ each iteration
% h0     = (Optional) Initial value of the solution




function out = atomic_LTI_BB(In,opt)
y = In.ym;
u = In.u;

tau = opt.tau;


if isempty(y)
    Ns = size(y, 1);
    ny = opt.ny;
    nu = opt.nu;
    F = {};
else
    [Ns,ny] = size(y);
    nu = size(u,2);

    %Toeplitz the input for fast multiplication
    F = cell(1,nu); %fft toeplitz of input, possibly trivial if time domain data is not used
    for ku = 1:nu    
       %Tu{ku} = toeplitz(u(:,ku) ,[u(1,ku);zeros(Ns-1,1)]);
       F{ku}  = toeplitzmultaux(u(:,ku),  [u(1,ku);zeros(Ns-1,1)]);
       %matrix-vector toeplitz multiplication:
       %y2 = toeplitzmult2(F, u);
    end

end


%Process the poles into groups for sum-of-norms regularization
p  = In.PoleArray;
ha = In.ImpRespArray;
f  = In.FreqRespArray;
W  = In.FreqWeight;
np = length(p);

g = In.PoleGroups;
g_rep_0 = permute(repmat(g', [1, nu, ny]), [3,2,1]);
g_rep = squeeze(reshape(g_rep_0, 1, [], 1));
g_hot = ind2vec(g); %one-hot encoding of groups
g_order = sum(g_hot, 2);
Ngroups = max(g);

%testing purposes only
%enforce all subsystems to have the same weight
%trying to figure out why reweighted heuristic increases error between
%iterations, if the setpoint is feasible on both atomic balls
SAME_WEIGHT = 0;

if isfield(In, 'PoleGroupWeights')
    gw =  In.PoleGroupWeights;
else
    if SAME_WEIGHT 
        gw = ones(size(g_order));
    else
        gw = g_order;	

        if p(1) == 1
            gw(1) = In.ConstWeight;
        end
    end
end

%deal with constant input

%g_offset = (0: (ny*nu-1))*np;
%g_offset = 0:(nu*ny-1);


w.groups = cell(Ngroups, 1);
%w.weights = zeros(Ngroups, 1);
w.weights = gw;
w.order = g_order;


for gi = 1:Ngroups    
    w.groups{gi} = find(g_rep  == gi);
end

%Formulate the least squares operators and paramters
%system to output wrt. input and its adjoint

b_time = reshape(y, Ns*ny, 1);

%Frequency response penalization
Wdim = length(size(W));
if isempty(W)
    Wdim = 0;
end


%check the b_freq calculations
operator = struct;
%if Wdim == 3
if ~isempty(W)
    if isvector(W)
        %same frequency weighting on all systems
        W = repmat(W, [1, np, nu]);
    end
    %IO (weighting function for each input/output pair)
    operator.A  = @(x) mimo_io_A(x, np, nu, ny, Ns, F, ha, f, W);
    operator.At = @(r) mimo_io_At(r,np, nu, ny, Ns, F, ha, f, W);
    %Y_rep = repmat(Y, 1, 1, nu);
    %Y_rep = W.*permute(Y_rep, [1, 3, 2]);
%     G_ref = zeros(size(W, 1),  nu, ny);
%     for i = 1:ny
%         for j = 1:nu
%             G_ref(:, j, i) = Y(:, i)./U(:, j);
%         end
%     end
    G_ref = W.*opt.FreqResponse;
    b_freq = complex_unfold(squeeze(reshape(G_ref, [], 1, 1)));
    %b_freq = complex_unfold(squeeze(reshape(Y_rep, [], 1, 1)));
    
    %b_freq = complex_unfold(kron(reshape(permute(Y, [2,1]), [], 1), ones(nu, 1)), 1);
% elseif Wdim == 2
%     %Output (weighting function for each output)
%     %don't use this
%     A  = @(x) mimo_output_A(x, np, nu, ny, Ns, F, ha, f, U, W);
%     At = @(r) mimo_output_At(r,np, nu, ny, Ns, F, ha, f, U, W);
%     %b_freq = complex_unfold(reshape(permute(Y, [2,1]), [], 1));
%     %b_freq = complex_unfold(reshape(W.*Y, [], 1));
%     b_freq = complex_unfold(reshape(W, [], 1));
else   
    %Time (no frequency penalization)
    operator.A  = @(x) mimo_A(x, np, nu, ny, Ns, F, ha);
    operator.At = @(r) mimo_At(r,np, nu, ny, Ns, F, ha);
    b_freq = [];
end

%b_freq = squeeze(reshape(G_ref, [], 1, 1));

%b_freq = repmat(Y, 1, 1, nu);
%

%I really hope this works
%want to copy Y, slicing by indices
%probably easier to write it out or keep everything as arrays
%or just let reshape take care of everything
%b_freq = kron(reshape(permute(Y, [2,1]), [], 1), [1;1;1]);
%b_freq = squeeze(reshape(repmat(permute(Y, [2,1]), 2,1,1), [], 1));
%
%b_freq = squeeze(reshape(b_freq, 1, 1 ,[] ));


b = [b_time; b_freq];

%Weighting function?



%something about reshaping?

BB_opt.num_var = nu*ny*np;
BB_opt.tau = tau;
BB_opt.w = w;
%BB_opt.delta = 0;
BB_opt.delta = 1e-3;
%BB_opt.delta = opt.delta;
%BB_opt.norm_type = 2;
%BB_opt.norm_type = Inf;
BB_opt.norm_type = In.NormType;
BB_opt.is_complex = 0;
% BB_opt.visualize = 0;
BB_opt.visualize_end = 0;
%BB_opt.DG_tol = 3e-3;
BB_opt.DG_tol = 1e-2;
BB_opt.ASQP = 0;

if isfield(In, 'warm_start')
     if isstruct(In.warm_start)
         BB_opt.warm_start = In.warm_start;
         BB_opt.warm_start.bash_manager.bagger.w = w;
     end     
     BB_opt.export_warm_start = 1;
end



tic
%Run the optimization routine

if opt.FCFW
    [x_final, S_final, c_final, run_log] = BB_operator(operator, b, BB_opt);
    out.Atoms = S_final;
    out.AtomCoeff = c_final;
else
    [x_final, run_log] = BB_forward(operator, b, BB_opt);
end


%[~, x_norm] = LMO_1d(x_final, opt.NormType, w);
%fprintf('Gap: %0.3e\n', opt.tau - x_norm)
out.time = toc;


out.Coeff0 = x_final;




out.run_log = run_log;


%output from data
Ax = operator.A(x_final);
%h_all = A(x_final)

%time
%y_time = b(1:(Ns*ny));
y_time = Ax(1:(Ns*ny));
out.y = reshape(y_time, Ns, ny);

%Frequency
%y_freq_real = b(Ns*ny + 1:end);
y_freq_real = Ax(Ns*ny + 1:end);
y_freq = complex_fold(y_freq_real, 1);

%if Wdim == 3
if ~isempty(W)
    out.f = reshape(y_freq, size(f, 1), nu, ny);
% elseif Wdim == 2
%     out.f = reshape(y_freq, Ns, ny);
else
    out.f = [];
end
%determine output
%x_coeff = reshape(full(x_final), np, nu, ny);
%x_coeff = reshape(full(x_final), ny, nu, np);

out.Coeff = cell(ny, nu);
out.h = cell(ny, nu);

for j = 1:nu
    for i = 1:ny
        %x_curr = sparse(squeeze(x_coeff(i, j, :)));
        ind_curr = sub2ind([ny, nu, np], i*ones(1, np), j*ones(1,np), 1:np);
        x_curr = x_final(ind_curr);
        out.Coeff{i, j} = x_curr;
        out.h{i, j} = ha*x_curr;
    end
end

out.iter = length(run_log.time);
%out.cost = norm(Ax - b)^2/2;
out.cost = run_log.error_list(end); %this is cheating, why are the values in disagreement?


%reweighting?
epsilon = 1e-4;
%delta = 0;
out.group_active =  [];
%out.PoleGroupWeights_old_all = gw;
%out.PoleGroupWeights_new= [];
weights_old = [];
weights_new = [];
weights_new_all = [];
x_max = [];
residues = [];

out.sys_out= 0;
out.sys_modes = {}; %In.sys_modes;
out.poles_active_ind = [];
out.poles_active = [];
z = tf('z', opt.SampleTime);
for gi = 1:Ngroups
    g_curr = w.groups{gi};    
    x_curr = x_final(g_curr);       
    x_max_curr = norm(x_curr, 'inf');        
    
    if x_max_curr
        x_max(end+1) = x_max_curr;
        weights_old(end+1) = gw(gi);
        if SAME_WEIGHT
            weights_new(end+1) = 1/(epsilon +  x_max_curr);
        else
            if gi == 1 && p(1) == 1
                weights_new(end+1) = In.ConstWeight/(epsilon +  x_max_curr);
            else    
                weights_new(end+1) = w.order(gi)/(epsilon +  x_max_curr);
            end
        end       
        poles_curr_ind = find(In.PoleGroups == gi);
        poles_curr = In.PoleArray(poles_curr_ind );
        
        out.group_active(end+1) = gi;
        out.poles_active_ind = [out.poles_active_ind poles_curr_ind];
        out.poles_active = [out.poles_active poles_curr];
        
        %Return the resultant system
        if opt.FormSystem
            x_curr_box = reshape(full(x_curr), ny, nu, []);
            %
            scales_curr = getScales(poles_curr, Ns);
            %use complex scales instead
            %scale_curr = getScales2(poles_curr(1), Ns);           


            if w.order(gi) == 1
                residue_curr = x_curr_box*scales_curr;
                sys_curr = residue_curr/(z - poles_curr);
            else

                res_cos =  x_curr_box(:, :, 1)*scales_curr(1);
                res_sin =  x_curr_box(:, :, 2)*scales_curr(2);
                      
                num_std  = res_cos - 1.0j*res_sin;
                num_conj = res_cos + 1.0j*res_sin;
                
                sys_std = num_std/(z - poles_curr(1));
                
                sys_conj = num_conj/(z - poles_curr(2));
                
                
                sys_curr = (sys_std + sys_conj);
                
            end
            
            sys_curr = zpk(sys_curr);
          
            out.sys_out = out.sys_out + sys_curr;
            out.sys_modes{end+1} = sys_curr;
        end

    end            
end

out.system_order =full(sum(w.order(out.group_active)));
out.tau = opt.tau;
weights_new_norm = weights_new * out.tau/ (weights_new*x_max');
out.PoleGroupWeights_new = weights_new_norm;
out.PoleGroupWeights_old = weights_old;
%weight_compare = [weights_old; weights_new_norm]*x_max';

out.norm_max = x_max;

% weights_new_all = ones(1, Ngroups)/delta;
% weights_new_all(out.group_active) = weights_new_norm;
% out.PoleGroupWeights_new_all = weights_new_all;

%normalize new set of weights (?)
%not sure if this is correct normalization
%out.PoleGroupWeights_new = out.PoleGroupWeights_new*...
%    length(out.PoleGroupWeights_new)/sum(out.PoleGroupWeights_new);

% out.poles_active_ind = any(In.PoleGroups == out.group_active', 1);
% out.poles_active = In.PoleArray(out.poles_active_ind)';   <-- THERES THE
% BS TRANSPOSE THAT IS SWITCHING MY SIGNS! (swear words go here)

out.w = w;
% new pole groups
[C, ia, ic] = unique(In.PoleGroups(out.poles_active_ind));
out.PoleGroups_new = ic';


out.y_orig = In.ym;
out.u = In.u;
%out.residues_active = residues;
    
%h_svd = svd(hankel_mo(h(2:min(1000,N),:)'));
%out.svd = h_svd(1:min(15,length(h_svd)))';
%assignin('base','out',out);

%out.error = norm(In.h-out.h)/norm(In.h);
end