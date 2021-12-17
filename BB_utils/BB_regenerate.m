function [x_rec, S_rec, c_rec ,run_log_rec] = BB_regenerate(x, opt, DG_tol)
%BB_REGENERATE If the existing set of atoms are blocking efficient progress
%(wasting atomic budget), pick a new set of atoms that equivalently
%describes the input data x. Part of the  'cheat to win' initiative
%
%Given x = c1' a1 with L1(c1) = tau, attempt to find a x = c2' a2 with
%L1(c2) <= tau. This regeneration arises from a bug somewhere in the group
%sparsity code, and should not be needed if the code routines works.
data = struct;
Id = @(xi) xi;
data.A = Id;
data.At = Id;
w = opt.w;
opt_rec = opt;
opt_rec.delta = 0;
opt_rec.regen_depth = 1;
if nargin < 3   
    opt_rec.DG_tol = 1e-6;
else
    opt_rec.DG_tol = DG_tol;
end
if isfield(opt_rec, 'warm_start')
    opt_rec = rmfield(opt_rec, 'warm_start');
end
%opt_rec.visualize = 1;
[x_rec, S_rec, c_rec ,run_log_rec] = BB_operator(data, x, opt_rec);
end

