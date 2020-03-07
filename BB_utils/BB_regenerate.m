function [x_rec, S_rec, c_rec ,run_log_rec] = BB_regenerate(x, opt)
%BB_REGENERATE If the existing set of atoms are blocking efficient progress
%(wasting atomic budget), pick a new set of atoms that equivalently
%describes the input data x. Part of the  'cheat to win' initiative
%
%

Id = @(xi) xi;
w = opt.w;
opt_rec = opt;
opt_rec.delta = 0;
if isfield(opt_rec, 'warm_start')
    opt_rec = rmfield(opt_rec, 'warm_start');
end
[x_rec, S_rec, c_rec ,run_log_rec] = BB_operator(Id, Id, x, opt_rec);
end

