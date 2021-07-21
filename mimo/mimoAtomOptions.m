classdef mimoAtomOptions <handle
   % System requirements and algorithm options for a MIMO LTI
   % identification.
   % Two stage algorithm: Randomize and Reweight
   properties        
       
       
      %% General MIMO Optimization 
      %group norm used in atomic penalty
      NormType(1,1)double{mustBePositive} = Inf;
      
      tau(1,1) double{mustBeNonnegative} = 1; % Atomic radius 
      delta(1,1)double{mustBeNonnegative} = 0; %Elastic net regularization
      
      %Frequency domain constraints and weights
      FreqWeight= [];     % Weights for frequency domain fitting
      FreqSample = [];    % frequencies to sample at
      FreqResponse = [];  % frequency response           
            
      %use Fully Corrective Frank Wolfe (1), or Forward Frank Wolfe (0)
      FCFW(1,1) logical = true;   
      
      %TODO: number of inputs and outputs (corrrect this)
      nu(1,1) double{mustBeNonnegative} = 1;
      ny(1,1) double{mustBeNonnegative} = 1;
       
            
      %% Randomize (add new poles)
       
      % Rational atom circle sector configured by radius and angle limits
      r1(1,1)double {mustBeNonnegative,mustBeLessThanOrEqual(r1,1)} = 0; % lower radius bound for rational atoms
      r2(1,1)double {mustBeNonnegative,mustBeLessThanOrEqual(r2,1)} = 1; % upper radius bound for rational atoms
      phi1(1,1)double {mustBeNonnegative,mustBeLessThanOrEqual(phi1,3.1415926535898)} = 0; % minimum angle, radians, for rational atoms
      phi2(1,1)double {mustBeNonnegative,mustBeLessThanOrEqual(phi2,6.2833)} = 2*pi; % maximum angle, radians, for rational atoms
      
      % Include +/- 1 as an atom, regardless of "Type" 
      IncludeConstant(1,1)logical = true; % include +/-1
      
      % Number of atoms
      NumAtoms(1,1)double{mustBePositive,mustBeInteger} = 1e3;
              
      %Number of rounds of randomization after the initial pick of poles
      RandomRounds(1,1)double{mustBeNonnegative,mustBeInteger}  = 0; 
            
      %% Reweight (sparsify existing set of poles)
      
      %Number of rounds of reweighted heuristic to sparsify the poles
      ReweightRounds(1,1)double{mustBeNonnegative,mustBeInteger} = 0;
      
      %tolerance to stop reweighting iterations
      ReweightTol(1,1)double{mustBeNonnegative} = 1e-5;       
            
      %% Output and Comparision 
      FormSystem(1,1)logical = false;     %Form a transfer function of the system      
      Compare(1,1)logical = true; %compare against benchmarks
   end
end
