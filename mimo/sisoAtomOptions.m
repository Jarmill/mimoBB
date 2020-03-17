classdef sisoAtomOptions <handle
   % System requirements and algorithm options for a SISO LTI
   % identification.
   properties
      tau = 1; % Atomic radius 
      
      % randomize atomic set over iterations
      Randomize(1,1)logical = false;
      
      %use Fully Corrective Frank Wolfe (1), or Forward Frank Wolfe (0)
      FCFW(1,1) logical = true;
      
      AwayStepOnNeedBasis

      % nature of atoms
      % rational: rational atomic set (traditional)
      % spline: first order spline atoms parameterized by 0<=alpha<=1
      % random: same as spline, but alpha is changed in each iteration
      Type(1,1)string{mustBeMember(Type,["rational","spline","random","TC"])} = "rational"; 

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
      Delta
      
      %Frequency domain parameters
      FreqWeight= [];     % Weights for frequency domain fitting
      FreqSample = [];    % frequencies to sample at
      FreqResponse = [];  % frequency response 
      
     
      % Trace weighting matrix when IsDet is FALSE. Specify as ny-by-ny psd
      % matrix (ny := number of model outputs)
      Wt = 1; 
      RandomRounds = 0; %Number of rounds of randomization after the initial pick of poles
      ReweightRounds = 0; %Number of rounds of reweighted heuristic to sparsify the poles
      ReweightTol = 1e-3; %tolerance to stop reweighting
      FormSystem = 0;     %Form a transfer function of the system
      % Pick whther the poles should be shared among all I/O pairs during
      % estimation iterations ("ss") or can be picked separately for each
      % pair ("tf"). 
      ModelType(1,1)string{mustBeMember(ModelType,["ss","tf"])} = "ss";
           
      %compare against benchmarks
      Compare(1,1)logical = true;
   end
end
