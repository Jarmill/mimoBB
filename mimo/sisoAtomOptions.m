classdef sisoAtomOptions <handle
   % System requirements and algorithm options for a SISO LTI
   % identification.
   properties
      tau = 1; % Atomic radius 
      AdaptiveTau(1,1) logical = false;
      
      % randomize atomic set over iterations
      Randomize(1,1)logical = false;
      SeedBasis(1,1)logical = false;
      
      AwayStepOnN
      % nature of atoms
      % rational: rational atomic set (traditional)
      % spline: first order spline atoms parameterized by 0<=alpha<=1
      % random: same as spline, but alpha is changed in each iteration
      Type(1,1)string{mustBeMember(Type,["rational","spline","random","TC"])} = "rational"; 
      
      % FOSS: Alpha is a scalar; MFOSS: Row vector
      % When Type="random", the specified Alpha value is used only for the
      % first iteration
      Alpha(1,:)double {mustBeNonnegative,mustBeLessThanOrEqual(Alpha,1)} = 0.6 % SS atom frequency constant 
      
      % Rational atom circle sector configured by radius and angle limits
      r1(1,1)double {mustBeNonnegative,mustBeLessThanOrEqual(r1,1)} = 0; % lower radius bound for rational atoms
      r2(1,1)double {mustBeNonnegative,mustBeLessThanOrEqual(r2,1)} = 1; % upper radius bound for rational atoms
      phi1(1,1)double {mustBeNonnegative,mustBeLessThanOrEqual(phi1,3.1415926535898)} = 0; % minimum angle, radians, for rational atoms
      phi2(1,1)double {mustBeNonnegative,mustBeLessThanOrEqual(phi2,6.2833)} = pi; % maximum angle, radians, for rational atoms
      
      % Include +/- 1 as an atom, regardless of "Type" 
      IncludeConstant(1,1)logical = true; % include +/-1
      
      % Number of atoms
      NumAtoms(1,1)double{mustBePositive,mustBeInteger} = 1e3;
      % Number of trace iterations
      MaxIterTrace(1,1)double{mustBePositive,mustBeInteger} = 200;
      % Number of det iterations
      MaxIterDet(1,1)double{mustBePositive,mustBeInteger} = 200;
      
      % Unit disk sampling. "uniform", or "normal"
      SamplingType(1,1)string{mustBeMember(SamplingType,["uniform","normal"])} = "uniform"; 
      SamplingFcn = []; % function handle to pick valid member atoms
      
      wn param.Continuous % Bounds on natural frequencies for a certain I/O pair 
      zeta param.Continuous % Bounds on damping 
      DCGain param.Continuous % Bound on dc gain
      MaxOvershoot param.Continuous % Bound on maximum overshoot
      Bandwidth param.Continuous % Bound on bandwidth
     
      %group norm used in atomic penalty
      NormType(1,1)double{mustBePositive} = Inf;
      
      %Frequency domain parameters
      FreqWeight= [];     % Weights for frequency domain fitting
      FreqSample = [];    % frequencies to sample at
      FreqResponse = [];  % frequency response 
      
      % Are we minimizing the determinant criterion det(EE') or the trace
      % criterion trace(EE')? Matters only in multi-output case. 
      IsDet = false; % trace or det criterion
      
      % Trace weighting matrix when IsDet is FALSE. Specify as ny-by-ny psd
      % matrix (ny := number of model outputs)
      Wt = 1; 
      RandomRounds = 0; %Number of rounds of randomization after the initial pick of poles
      ReweightRounds = 0; %Number of rounds of reweighted heuristic to sparsify the poles
      FormSystem = 0;     %Form a transfer function of the system
      % Pick whther the poles should be shared among all I/O pairs during
      % estimation iterations ("ss") or can be picked separately for each
      % pair ("tf"). 
      ModelType(1,1)string{mustBeMember(ModelType,["ss","tf"])} = "ss";
      
      % FW scheme: 
      %  F:    Forward-only; one subsystem at a time
      %  FA:   Forward and Away; one subsystem at a time
      %  FC:   Fully Corrective; one subsystem at a time
      %  JD:   Joint Descent; Forward selection only; all subsystems
      %        together
      %  JDA:  Joint Descent; Forward and Away steps; all subsystems
      %        together
      % none:  Do not use FW.
      %  
      FW(1,1) string {mustBeMember(FW,["F","FA","FC","JD","JDA","none"])} = "FA";
      
      % Descent direction determination
      % gn: Gauss-Newton
      % gna: Adaptive Gauss-Newton
      % grad: Steepest Descent
      SearchMethod(1,1)string{mustBeMember(SearchMethod,["gn","grad","gna","lm"])} = "grad";
      
      % Show pole-zero map during estimation (slows things down)
      ShowProgressPlot(1,1)logical = true;
      
      %compare against benchmarks
      Compare(1,1)logical = true;
      
      % lsqlin options
      LsqlinOpt = optimoptions('lsqlin','Display','none');
      
      % Prefilter
      PF = 1;
   end
end
