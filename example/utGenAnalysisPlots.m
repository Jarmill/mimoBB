function utGenAnalysisPlots(out,sys)
% Quality analysis

z = out.Results.z;
syse = out.Results.syse; % original estimated (FIR)
%sysr = ss(out.Results.sysr); % reduced by balred (SS)
sysr2 = ss(out.Results.sysr2); % reduced by n4sid
sys_tfest_fd = out.Results.sys_tfest_fd;
sys_tfest_td = out.Results.sys_tfest_td;
sys_ssest_fd = out.Results.sys_ssest_fd;
sys_ssest_td = out.Results.sys_ssest_td;
sys_ssregest_td = out.Results.sys_ssregest_td;

[~,T] = impulse(sys);


% Pole location comparison
fprintf('Pole locations (true, estimated)\n')
[sort(pole(sys)), sort(pole(sysr2))]

fprintf('model output vs. measured output comparison\n')
fig1 = figure;
compare(z,sys,sys_ssest_td,sys_ssregest_td,sysr2,'init','z')
%compare(z,sys_tfest_fd,sys_tfest_td,sys_ssest_fd,sys_ssest_td,sys_ssregest_td,syse,sysr,'init','z')
title('Response fit comparison (100% is best)')
shg

fprintf('impulse response comparisons\n')
fig2 = figure;
impulse(sys,sys_ssest_td,sys_ssregest_td,sysr2,T(end))
title('My best result (sysr) impulse response compared to the true one')
shg
fig3 = figure;
%impulse(sys,sys_tfest_fd,sys_tfest_td,sys_ssest_fd,sys_ssest_td,sys_ssregest_td,syse,sysr2,T(end))
impulse(sys,sys_ssest_td,sys_ssregest_td,syse,sysr2,T(end))

title('All impulse responses compared to the true one')
shg
fig4 = figure;
%bodemag(sys,sys_tfest_fd,sys_tfest_td,sys_ssest_fd,sys_ssest_td,sys_ssregest_td,syse,sysr2)
bodemag(sys,sys_ssest_td,sys_ssregest_td,syse,sysr2)

title('Bode mag')
shg
