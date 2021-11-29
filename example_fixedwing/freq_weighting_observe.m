%identification of modes of a flexible wing aircraft system

%1 input:       driving spring at center of mass
%20 outputs:    accelerometers on wing

%10 experiments. Each experiment records 2 accelerometer time traces.
%Inputs are a frequency sweep from 3-35 Hz.

%https://www.mathworks.com/help/ident/ug/modal-analysis-of-a-flexible-flying-wing-aircraft.html

SOLVE = 1;
PLOT = 1;

load('flexwing_freq.mat');

rng(77, 'twister')

ny = 20;
nu = 1;

%do not use time samples yet.
Ns = 0;

%only frequency samples
Nf = length(Gs.Frequency);

%% pose sysid problem

%weighting matrix
%change to the square root of trace 11, like the demo later
W = ones(ny, nu, Nf);

opt = mimoAtomOptions;
opt.FreqWeight = W;
opt.FreqResponse = Gs.ResponseData;
opt.FreqSample = Gs.Frequency;
opt.Compare = 0;
opt.tau = 2;
opt.RandomRounds = 3;
opt.ReweightRounds = 10;
opt.NumAtoms = 200;
opt.NormType = Inf;
opt.FormSystem = 1;
opt.FCFW = 1; 
opt.SampleTime = 1/2000;

opt.IncludeConstant = 0; %disable later

%sector bounds
opt.r1 = 0.4; %making this up
opt.phi2 = 2*pi*(120)*Gs.Ts;

opt.ny = ny;
opt.nu = nu;

if SOLVE
    %no time domain, only frequency domain
[out, out_random] = atomic_LTI_iteration([], [], opt);

sys_all = 0;
for i = 1:length(out.sys_modes)
    sys_all = sys_all + out.sys_modes{i};
end

end