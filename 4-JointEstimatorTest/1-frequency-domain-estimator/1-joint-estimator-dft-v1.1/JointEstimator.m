function [xBest, yBest, info] = JointEstimator(xn, Fs, paramRange, ...
    options, optionsPso, optionsGrad)

%
% Joint estimator of frequency and phase of sinusoid
% Frequency domain, DFT-based estimation
% Correlation-based method
% Using self-coded PSO algorithms
% Using MATLAB fmincon as gradient method
% 
% Input arguments:
%   @xn         : Signal to be estimated
%   @Fs         : Sampling rate
%   @paramRange : Estimation parameter range [fLb, fUb, pLb, pUb]
%   @options    : Options for the whole estimator
%   @optionsPso : Options for PSO algorithm
%   @optionsGrad: Options for gradient optimization algorithm
%
% Output arguments:
%   @xBest  : Optimal point (variable)
%   @fBest  : Optimal value of object function
%   @info   : Information of the optimization process
%
% Author        : Zhiyu Shen @Nanjing University
% Establish Date: Aug 3, 2022
% Revised Data  : Dec 9, 2022
%


%% Step 1: Parameter Configuartion

% Input Vector Size Validation
n = size(xn, 1);
if n ~= 1
    error('Input signal must be in a row vector!');
end

% Defult set of options for the whole estimator
default.maxIter = 5;            % Maximum iteration times
default.display = 3;            % Print iteration progress out on the screen

% Set options according to user inputs
options = MergeOptions(default, options);

% Assign some paramters
maxIter = options.maxIter;

% Display options
% 0: Display each iteration in joint estimator
% 1: Display each iteration in particle swarm optimization
% 2: Display each iteration in conjugate gradient algorithm
if options.display == 1
    optionsPso.display = 'iter';
    optionsGrad.Display = 'off';
elseif options.display == 2
    optionsPso.display = 'none';
    optionsGrad.Display = 'iter';
elseif options.display == 3
    optionsPso.display = 'none';
    optionsGrad.Display = 'off';
end


%% Step 2: Preliminary Estimation (DTFT Peak Search)

% Add window and perform DFT on the test signal
nSam = length(xn);                              % Number of original signal samples
nFFT = 2^nextpow2(nSam);                        % Number of FFT points
winSig = ones(1,nFFT);                          % Window signal (Rectangular Window)
% idxWin = 0 : 1 : nFFT-1;
% winSig = 0.54 - 0.46*cos(2*pi*idxWin/nFFT);     % Window signal (Hamming Window)
xnWin = [xn, zeros(1,nFFT-nSam)].*winSig;       % Zero padding and add window
xnFFT = fft(xnWin,nFFT);

% Compute the frequency spectrum of test signal
Xn1 = abs(xnFFT/nFFT);
Xw = Xn1(1:nFFT/2);
Xw(2:end-1) = 2*Xw(2:end-1);

% Compute mean and variance of test signal
miu0 = sum(Xw);
sigma0 = sqrt(sum((Xw-miu0).^2)/nFFT);

% Compute signal information for correlation computation
Ct = (Xw-miu0) ./ sigma0;
Ct = [Ct, nSam, nFFT];

% Estimate with DFT peak search
xPeak = RoughEstimate(xnFFT, Xw, nFFT, Fs, Ct);
fPeak = xPeak(1);
pPeak = xPeak(2);

%%% Initialization

% Allocate memory for info
info.globalBestFreq = zeros(1, maxIter);        % Global best frequency value of current iteration
info.globalBestPha  = zeros(1, maxIter);        % Global best phase value of current iteration
info.globalBestFval = zeros(1, maxIter);        % Global optimal objective function value of current iteration
info.bestFreq       = zeros(1, maxIter);        % Final best frequency value of current iteration
info.bestPha        = zeros(1, maxIter);        % Final best phase value of current iteration
info.bestFval       = zeros(1, maxIter);        % Optimal objective function value of current iteration
info.bestFreqGrad   = zeros(1, maxIter);        % Frequency component of gradient of current iteration
info.bestPhaGrad    = zeros(1, maxIter);        % Phase component of gradient of current iteration
info.iterationTime  = zeros(1, maxIter);        % Time spend on each iteration
info.meanTime       = [];                       % Mean time spend on each iteration
info.iteration      = 1 : maxIter;

% Setup display header
if options.display == 0
    fprintf('\n                                                     frequency      phase\n');
    fprintf(  'Iteration      frequency      phase        f(x)      gradient      gradient\n');
end


%%% Search process

% Search range
fLb = fPeak*(1-0.2);
fUb = fPeak*(1+0.2);
pLb = pPeak*(1-0.2);
pUb = pPeak*(1+0.2);
xLb = [fLb, pLb];
xUb = [fUb, pUb];

% Iteration
xBest = zeros(1, 2);
yBest = 8;
for iter = 1 : maxIter

    startTime = tic;
    
    % Global search with random start
    nvars = 2;
    [xGlobal, yGlobal, ~] = ParticleSwarmOptim(Ct, Fs, ...
        nvars, xLb, xUb, optionsPso);

    % Local search
    fLbLoc = max(fLb, xGlobal(1)-0.05);
    fUbLoc = min(fUb, xGlobal(1)+0.05);
    pLbLoc = xGlobal(2)-pi/50;
    pUbLoc = xGlobal(2)+pi/50;
    lb = [fLbLoc, pLbLoc];
    ub = [fUbLoc, pUbLoc];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];
    fun = @(X)ObjFunFreq(X, Ct, Fs);
    [xIter, yIter, ~, ~, ~, gIter] = fmincon(fun, xGlobal, A, b, Aeq, ...
        beq, lb, ub, nonlcon, optionsGrad);
    

    % Log Data
    info.globalBestFreq = xGlobal(1);
    info.globalBestPha  = xGlobal(2);
    info.globalBestFval = yGlobal;
    info.bestFreq(iter) = xIter(1);
    info.bestPha(iter) = xIter(2);
    info.bestFval(iter) = yIter;
    info.bestFreqGrad(iter) = gIter(1);
    info.bestPhaGrad(iter) = gIter(2);
    info.iterationTime(iter) = toc(startTime);
    
    % Whether new iteration is better
    if yIter < yBest
        xBest = xIter;
        yBest = yIter;
    end % end: if
    
    % Print search result of current iteration
    if options.display == 0
        fprintf('%5.0d          %.3f Hz    %.3f rad       %.3f      %.3f         %.3f\n', ...
        iter, xIter(1), xIter(2), yIter, ...
        abs(gIter(1)), abs(gIter(2)));
    end % end: if

end % end: for

info.meanTime = sum(info.iterationTime) / maxIter;

end % end: function JointEstimator



%%%% Function "RoughEstimate"

function xBest = RoughEstimate(xnFFT, Xw, nFFT, Fs, Ct)
%
% Frequency estimator with DFT peak search
% Dichotomous method
% 
% Input arguments:
%   @xnFFT : DFT value of test signal
%   @Xw    : Frequency spectrum value of test signal
%   @nFFT  : Number of FFT points
%   @Fs    : Sampling frequency
%   @Ct    : Covariance information of test signal's frequency spectrum
%
% Output arguments:
%   @xPeak: Frequency and phase result of the DFT peak search
%

% Find the peak of DFT
[~, idxPeak] = max(Xw);
p1 = angle(xnFFT(idxPeak));

% Compute the frequency value of the peak and the secondary peak
f1 = (idxPeak-2)*Fs/nFFT;
f2 = idxPeak*Fs/nFFT;

% Compute the correlation coefficient at 
% the left and right end of the frequency search range
r1 = ObjFunFreq([f1 p1], Ct, Fs); 
r2 = ObjFunFreq([f2 p1], Ct, Fs); 

% Find the maximum correlation coefficient with dichotomous method
while f2-f1 > 1e-9
    fm = (f1+f2)/2;
    rm = ObjFunFreq([fm p1], Ct, Fs); 
    if r2 < r1 
        f1 = fm;
        r1 = rm;
    else
        f2 = fm;
        r2 = rm;
    end
end
fPeak = (f1+f2)/2;

% Estimate phase roughly
if idxPeak == 1
    f1 = (idxPeak-1)*Fs/nFFT;
    f2 = idxPeak*Fs/nFFT;
    idx1 = idxPeak;
    idx2 = idxPeak+1;
else
    if Xw(idxPeak-1) > Xw(idxPeak+1)
        f1 = (idxPeak-2)*Fs/nFFT;
        f2 = (idxPeak-1)*Fs/nFFT;
        idx1 = idxPeak-1;
        idx2 = idxPeak;
    else
        f1 = (idxPeak-1)*Fs/nFFT;
        f2 = idxPeak*Fs/nFFT;
        idx1 = idxPeak;
        idx2 = idxPeak+1;
    end
end

% Compute the angle value at peak and correct it
p1 = angle(xnFFT(idx1));
p2 = angle(xnFFT(idx2));
if abs(p1-p2) > 5
    if p1 > p2
        p2 = p2 + 2*pi;
        p0 = angle(xnFFT(idx2)) + (fPeak-f2)*(p1-p2)/(f1-f2);
    else
        p1 = p1 + 2*pi;
        p0 = angle(xnFFT(idx2)) + (fPeak-f2)*(p1-p2)/(f1-f2);
    end
else
    p0 = p2 + (fPeak-f2)*(p1-p2)/(f1-f2);
end
pPeak = mod(p0, 2*pi);

% Generate output vector
xBest = [fPeak, pPeak];

end % end: function RoughEstimate



%%%% Function "ParticleSwarmOptim"

function [xBest, fBest, info] = ParticleSwarmOptim(Ct, ...
    Fs, nvars, xLb, xUb, options)
%
% Intelligent optimization algorithm called Particle Swarm Optimization
% Adopted in global search to get a rough estimation of the optimal solution
%
% Input arguments:
%   @Ct     : Covariance information of sequence to be estimated
%   @Fs     : Sampling rate
%   @nVar   : Variable dimension
%   @xLb    : Lower bound of variables (1*nvars)
%   @xUb    : Upper bound of variables (1*nvars)
%   @options: Optimization options, for more details see 'Option Defult
%             Set' in 'Preparation' part
%
% Output arguments:
%   @xBest         : Optimal point (variable)
%   @fBest         : Optimal value of object function
%   @totalTime     : Total time for optimization
%   @totalIteration: Total iteration times
%

%%% Preparation

% Input vector size validation
% xLb, xUb: nV*1 matrix
[n, m] = size(xLb);
if (n ~= 1) || (m ~= nvars)
    error('xLb is not of a valid size!')
end % end: if
[n, m] = size(xUb);
if (n ~= 1) || (m ~= nvars)
    error('xUb is not of a valid size!')
end % end: if

% Set options according to user inputs
if nargin == 6
    userOptions = options;
else
    userOptions = [];
end % end if
options = SetOptions(userOptions, nvars, Ct, Fs);

% Assign parameters
numParticles = options.swarmSize;                   % Particle(population) size
cSelf = options.selfAdjustment;                     % Self adjustment weight
cSocial = options.socialAdjustment;                 % Social adjustment weight
minNeighborhoodSize = max(2, ...
    floor(numParticles*options.minNeighborFrac));   % Minimum neighborhood size
minInertia = options.inertiaRange(1);               % Minimum inertia
maxInertia = options.inertiaRange(2);               % Maximum inertia
xLbMat = repmat(xLb, numParticles, 1);              % Particle position lower bound matrix
xUbMat = repmat(xUb, numParticles, 1);              % Particle position upper bound matrix


%%% Initialization

% Create initial state: particle positions & velocities, fvals, status data
state = MakeState(nvars, xLbMat, xUbMat, options);
bestFval = min(state.fvals);
% Create a vector to store the last StallIterLimit bestFvals.
% bestFvalsWindow is a circular buffer, so that the value from the i'th
% iteration is stored in element with index mod(i-1,StallIterLimit)+1.
bestFvalsWindow = nan(options.stallIterLimit, 1);

% Initialize adaptive parameters:
%   Initial inertia = maximum * magnitude * inertia
%   Initial neighborhood size = minimum neighborhood size
adaptiveInertiaCounter = 0;
if all(options.inertiaRange >= 0)
    adaptiveInertia = maxInertia;
elseif all(options.inertiaRange <= 0)
    adaptiveInertia = minInertia;
end % end: if
adaptiveNeighborhoodSize = minNeighborhoodSize;

% Setup display header
if options.verbosity > 1
    fprintf('\n                                 Best            Mean     Stall\n');
    fprintf(  'Iteration     f-count            f(x)            f(x)    Iterations\n');
    fprintf('%5.0f         %7.0f    %12.4g    %12.4g    %5.0f\n', ...
        0, state.funEval, bestFval, mean(state.fvals), 0);
end


%%% Iteration

pIdx = 1 : numParticles;                % Particle index
exitFlag = [];
% Run the main loop until some exit condition becomes true
while isempty(exitFlag)

    state.iteration = state.iteration + 1;

    bestNeighborIdx = GenerateBestNeighborIdx(state, ...
        adaptiveNeighborhoodSize, numParticles);

    % Update the velocities
    state.velocities(pIdx,:) = UpdateVelocities(state, adaptiveInertia, ...
        bestNeighborIdx, cSelf, cSocial, pIdx, nvars);

    % Update the positions
    [state.positions(pIdx, :), tfInvalid] = UpdatePositions(state, ...
        xLbMat, xUbMat, pIdx, numParticles, nvars);

    % For any particle on the boundary, enforce velocity = 0.
    if any(tfInvalid(:))
        state.velocities(tfInvalid) = 0;
    end
    
    % Update the objective function values
    state.fvals = ObjFunFreq(state.positions, Ct, Fs);

    % Update state with best fvals and best individual positions
    state = UpdateState(state, numParticles, pIdx);
    
    bestFvalsWindow(1 + mod(state.iteration-1, ...
        options.stallIterLimit)) = min(state.individualBestFvals);
    
    % Update inertia factor
    [state, adaptiveInertiaCounter, bestFval, adaptiveNeighborhoodSize, ...
        adaptiveInertia] = UpdateInertia(state, minInertia, ...
        maxInertia, bestFval, adaptiveInertiaCounter, ...
        adaptiveNeighborhoodSize, adaptiveInertia, numParticles, ...
        minNeighborhoodSize);

    % check to see if any stopping criteria have been met
%     [exitFlag, stopCondition] = StopParticleswarm(options, state, ...
%         bestFvalsWindow);

    [exitFlag, ~] = StopParticleswarm(options, state, ...
        bestFvalsWindow);

end % End while loop

% Find and return the best solution
[fBest, indexBestFval] = min(state.individualBestFvals);
xBest = state.individualBestPositions(indexBestFval,:);

% Generate output information
info.totalTime = toc(state.startTime);
info.totalIteration = state.iteration;

end % end: function ParticleSwarmOptim



%%%% Function "SetOptions"

function options = SetOptions(userOptions, nvars, Ct, Fs)
%
% Set optimization options
%
% Input arguments:
%   @userOptions: Options defined by user
%   @nvars      : Number of variables
%   @Ct         : Covariance information of sequence to be estimated
%   @Fs         : Sampling rate
%
% Output arguments:
%   @options: Options for optimization
%

% Option Defult Set
default.inertiaRange       = [0.1, 1.1];            % Inertia range
default.selfAdjustment     = 1.49;                  % Self adjustment factor
default.socialAdjustment   = 1.49;                  % Social adjustment factor
default.swarmSize          = min(100, 10*nvars);    % Number of particles
default.initialSwarm       = [];                    % Initial particle swarm
default.initialSwarmSpan   = 2000;                  % Initial swarm span
default.maxGeneration      = min(100, 200*nvars);   % Maximum number of generatons
default.maxTime            = inf;                   % Maximum time the algorithm runs
default.minNeighborFrac    = 0.25;                  % Minimum neighborhood size fraction
default.objectiveLimit     = -inf;                  % Minimum objective function value desired
default.tolFunValue        = 1e-9;                  % Termination tolerance on function value
default.stallIterLimit     = 20;                    % Maximum number of stalled iterations
default.stallTimeLimit     = inf;                   % Maximum time of stalled iterations
default.display            = 'none';                % Whether iteration progress printed on the screen
default.displayInterval    = 1;                     % Iteration interval printed on the screen
default.displaySectionSize = 20;                    % Number of iterations displayed in a section
default.covarianceMatrix   = Ct;                    % Covariance information of sequence to be estimated
default.samplingFrequency  = Fs;                    % Sampling rate

% Merge user defined options with default ones
if ~isempty(userOptions)
    options = MergeOptions(default, userOptions);
else
    options = default;
end % end: if

% Determine the verbosity
switch  options.display
    case {'off','none'}
        options.verbosity = 0;
    case 'final'
        options.verbosity = 1;
    case 'iter'
        options.verbosity = 2;
end % end: switch

end % end: function SetOptions



%%%% Function "MakeState"

function state = MakeState(nvars, xLbMat, xUbMat, options)
%
% Create an initial set of particles and objective function values
%
% Input arguments:
%   @nvars  : Number of variables
%   @xLb    : Lower bound of particle position
%   @xUb    : Upper bound of particle position
%   @options: Optimization options
%
% Output arguments:
%   @state: State struct
%

% makeState needs the vector of bounds, not the expanded matrix.
xLb = xLbMat(1, :);
xUb = xUbMat(1, :);

% A variety of data used in various places
state = struct;
state.iteration = 0;                % Current generation counter
state.startTime = tic;              % Tic identifier
state.stopFlag = false;             % OutputFcns flag to end the optimization
state.lastImprovement = 1;          % Generation stall counter
state.lastImprovementTime = 0;      % Stall time counter
state.funEval = 0;                  % Number of objective function evaluations
numParticles = options.swarmSize;    % Number of particles

% If InitialSwarm is partly empty use the creation function to generate
% population (CreationFcn can utilize InitialSwarm)
if numParticles ~= size(options.initialSwarm,1)
    swarmStruct = struct;
    swarmStruct.xLb = xLb;
    swarmStruct.xUb = xUb;
    swarmStruct.nvars = nvars;
    swarmStruct.options = options;
    state.positions = SwarmCreationUniform(swarmStruct);
else
    state.positions = options.initialSwarm;
end

% Enforce bounds
if any(any(state.positions < xLbMat)) || any(any(state.positions > xUbMat))
    state.positions = max(xLbMat, state.positions);
    state.positions = min(xUbMat, state.positions);
end

% Initialize velocities by randomly sampling over the smaller of initSwarmSpan or ub-lb
% min will be initSwarmSpan if either lb or ub is not finite
vMax = min(xUb-xLb, options.initialSwarmSpan);
state.velocities = repmat(-vMax, numParticles, 1) + ...
    repmat(2*vMax, numParticles, 1) .* rand(numParticles, nvars);

% Calculate the objective function value for all particles.
Ct = options.covarianceMatrix;
Fs = options.samplingFrequency;
fvals = ObjFunFreq(state.positions, Ct, Fs);
state.fvals = fvals;
state.funEval = numParticles;

state.individualBestFvals = state.fvals;
state.individualBestPositions = state.positions;

end % end: function MakeState



%%%% Function "SwarmCreationUniform"

function swarm = SwarmCreationUniform(swarmStruct)
%
% Creates the initial positions for PSO algorithm
%
% Input arguments:
%   @swarmStruct: State struct for swarm creation
%
% Output arguments:
%   @swarm: Initial particle swarm (numParticle*nvars)
%

% Assign parameters according to swarm creation struct
nvars = swarmStruct.nvars;
options = swarmStruct.options;
xLb = swarmStruct.xLb;
xUb = swarmStruct.xUb;

% Assign parameters according to options struct
numParticles = options.swarmSize;
numInitPositions = size(options.initialSwarm, 1);
numPositionsToCreate = numParticles - numInitPositions;

% Initialize particles to be created
swarm = zeros(numParticles, nvars);

% Use initial particles provided already
if numInitPositions > 0
    swarm(1:numInitPositions, :) = options.initialSwarm;
end % end: if

% Create remaining particles, randomly sampling within lb and ub
span = xUb - xLb;
swarm(numInitPositions+1:end, :) = repmat(xLb, numPositionsToCreate, 1) + ...
    repmat(span, numPositionsToCreate, 1) .* rand(numPositionsToCreate, nvars);

end % end: function SwarmCreationUniform



%%%% Function "GenerateBestNeighborIdx"

function bestNeighborIdx = GenerateBestNeighborIdx(state, ...
    adaptiveNeighborhoodSize, numParticles)
%
% Generate best neighborhood index
% The best particle in random neighborhood
% The size is controlled by the adaptiveNeighborhoodSize parameter
%
% Input arguments:
%   @state                   : State struct of optimization
%   @adaptiveNeighborhoodSize: Number of neighbors
%   @numParticles            : Number of particles
%
% Output arguments:
%   @bestNeighborIdx: Index of neighbor with best object function value
%

neighborIdx = zeros(numParticles, adaptiveNeighborhoodSize);
neighborIdx(:, 1) = 1 : numParticles;       % First neighbor is self
for i = 1 : numParticles
    % Determine random neighbors that exclude the particle itself,
    % which is (numParticles-1) particles
    neighbors = randperm(numParticles-1, adaptiveNeighborhoodSize-1);
    % Add 1 to indicies that are >= current particle index
    iShift = neighbors >= i;
    neighbors(iShift) = neighbors(iShift) + 1;
    neighborIdx(i, 2:end) = neighbors;
end % end: for

% Identify the best neighbor
[~, bestRowIndex] = min(state.individualBestFvals(neighborIdx), [], 2);

% Create the linear index into neighborIndex
bestLinearIdx = (bestRowIndex.'-1).*numParticles + (1:numParticles);
bestNeighborIdx = neighborIdx(bestLinearIdx);

end % end: function GenerateBestNeighborIdx



%%%% Function "UpdateVelocities"

function newVelocities = UpdateVelocities(state, adaptiveInertia, ...
    bestNeighborIndex,cSelf,cSocial,pIdx,nvars)
%
% Update the velocities of particles with indices pIdx
%
% Input arguments:
%   @state            : State struct of optimization
%   @adaptiveInertia  : Current inertia value
%   @bestNeighborIndex: Current best neighbor's index for each particle
%   @cSelf            : Current self adjustment factor
%   @cSocial          : Current social adjustment factor
%   @pIdx             : Particle index
%   @nvars            : Number of variables
%
% Output arguments:
%   @newVelocities: Updated particle velocities
%

% Generate random number distributions for self and social components
randSelf = rand(numel(pIdx), nvars);
randSocial = rand(numel(pIdx), nvars);

% Fetch old velocities from state structure
oldVelocities = state.velocities(pIdx, :);

% Update rule
newVelocities = adaptiveInertia*oldVelocities + ...
    cSelf*randSelf.*(state.individualBestPositions(pIdx,:)-state.positions(pIdx,:)) + ...
    cSocial*randSocial.*(state.individualBestPositions(bestNeighborIndex(pIdx), :)-state.positions(pIdx,:));

% Find infinite velocities, replace them with old velocities
tfInvalid = ~all(isfinite(newVelocities), 2);
newVelocities(tfInvalid) = oldVelocities(tfInvalid);

end % end: function UpdateVelocities



%%%% Function "UpdatePositions"

function [newPositions, tfInvalid] = UpdatePositions(state, xLbMat, ...
    xUbMat, pIdx, numParticles, nvars)
% 
% Update positions of particles with indices pIdx.
% 
% Input arguments:
%   @state       : State struct of optimization
%   @xLb         : Lower bound of particle position
%   @xUb         : Upper bound of particle position
%   @pIdx        : Particle index
%   @numParticles: Number of particles
%   @nvars       : Number of variables
%
% Output arguments:
%   @newPositions: Updated particle position
%   @tfInvalid   : Indicator flag for whether the particle positions 
%                  exceed bounds
%

newPositions = state.positions(pIdx, :) + state.velocities(pIdx, :);

% Remove positions if infinite.
tfInvalid = any(~isfinite(newPositions), 2);
tfInvalidFull = false(numParticles, 1);
tfInvalidFull(pIdx) = tfInvalid;
newPositions(tfInvalid, :) = state.positions(tfInvalidFull, :);

% Enforce bounds on positions and return logical array to update velocities
% where position exceeds bounds.
tfInvalidLb = newPositions < xLbMat(pIdx,:);
if any(tfInvalidLb(:))
    tfInvalidLBFull = false(numParticles, nvars);
    tfInvalidLBFull(pIdx, :) = tfInvalidLb;
    newPositions(tfInvalidLb) = xLbMat(tfInvalidLBFull);
    tfInvalid = tfInvalidLBFull;
else
    tfInvalid = false(numParticles,nvars);
end
tfInvalidUb = newPositions > xUbMat(pIdx,:);
if any(tfInvalidUb(:))
    tfInvalidUBFull = false(numParticles, nvars);
    tfInvalidUBFull(pIdx, :) = tfInvalidUb;
    newPositions(tfInvalidUb) = xUbMat(tfInvalidUBFull);
    tfInvalid = tfInvalid | tfInvalidUBFull;
end

end % end: function UpdatePositions



%%%% Function "UpdateState"

function state = UpdateState(state, numParticles, pIdx)
% 
% Update best fvals and best individual positions
%
% Input arguments:
%   @state       : State struct
%   @numParticles: Number of particles
%   @pIdx        : Particles index
%
% Output arguments:
%   @state: State struct
%

state.funEval = state.funEval + numel(pIdx);

% Remember the best fvals and positions for this block
tfImproved = false(numParticles, 1);
tfImproved(pIdx) = state.fvals(pIdx) < state.individualBestFvals(pIdx);
state.individualBestFvals(tfImproved) = state.fvals(tfImproved);
state.individualBestPositions(tfImproved, :) = state.positions(tfImproved, :);

end % end: function UpdateState



%%%% Function "UpdateInertia"

function [state, adaptiveInertiaCounter, bestFval, ...
    adaptiveNeighborhoodSize, adaptiveInertia] = UpdateInertia(state, ...
    minInertia, maxInertia, bestFval, adaptiveInertiaCounter, ...
    adaptiveNeighborhoodSize, adaptiveInertia, numParticles, minNeighborhoodSize)
% 
% Keep track of improvement in bestFvals and update the adaptive
% parameters according to the approach described in S. Iadevaia et
% al. Cancer Res 2010;70:6704-6714 and M. Liu, D. Shin, and H. I.
% Kang. International Conference on Information, Communications and
% Signal Processing 2009:1-5.
%
% Input arguments:
%   @state                   : State struct
%   @minInertia              : Minimum inertia factor
%   @maxInertia              : Maximum inertia factor
%   @bestFval                : Current best objective function value
%   @adaptiveInertiaCounter  : Current inertia stallation counter
%   @adaptiveNeighborhoodSize: Current number of neighbors
%   @adaptiveInertia         : Current inertia factor
%   @numParticles            : Number of particles
%   @minNeighborhoodSize     : Minumum number of neighbors
%
% Output arguments:
%   @state                   : State struct
%   @adaptiveInertiaCounter  : Updated inertia stallation counter
%   @bestFval                : Updated best objective function values
%   @adaptiveNeighborhoodSize: Updated number of neighbors
%   @adaptiveInertia         : Updated inertia factor
%

% Update global best particle
newBest = min(state.individualBestFvals);
if isfinite(newBest) && newBest < bestFval
    bestFval = newBest;
    state.lastImprovement = state.iteration;
    state.lastImprovementTime = toc(state.startTime);
    adaptiveInertiaCounter = max(0, adaptiveInertiaCounter-1);
    adaptiveNeighborhoodSize = minNeighborhoodSize;
else
    adaptiveInertiaCounter = adaptiveInertiaCounter+1;
    adaptiveNeighborhoodSize = min(numParticles, ...
        adaptiveNeighborhoodSize + minNeighborhoodSize);
end % end: if

% Update the inertia coefficient, enforcing limits (Since inertia
% can be negative, enforcing both upper *and* lower bounds after
% multiplying.)
if adaptiveInertiaCounter < 2
    adaptiveInertia = max(minInertia, min(maxInertia, 2*adaptiveInertia));
elseif adaptiveInertiaCounter > 5
    adaptiveInertia = max(minInertia, min(maxInertia, 0.5*adaptiveInertia));
end % end: if

end % end: function UpdateInertia



%%%% Function "StopParticleswarm"

function [exitFlag, reasonToStop] = StopParticleswarm(options, state, ...
    bestFvalsWindow)
% 
% Function handling conditions when iteration stops
% 
% Input arguments:
%   @options        : Optimization options
%   @state          : State struct
%   @bestFvalsWindow: Circular buffer storing previous iteration value
% 
% Output arguments:
%   @exitFlag    : Exit flag passed to outer loop
%   @reasonToStop: Iteration stop condition
% 

iteration = state.iteration;

iterationIdx = 1 + mod(iteration-1, options.stallIterLimit);
bestFval = bestFvalsWindow(iterationIdx);
if options.verbosity > 1 && ...
        mod(iteration, options.displayInterval)==0 && ...
        iteration > 0
    funEval  = state.funEval;
    meanFval = mean(state.fvals);
    stallGen = iteration - state.lastImprovement;
    fprintf('%5.0f         %7.0f    %12.4g    %12.4g    %5.0f\n', ...
        iteration, funEval, bestFval, meanFval, stallGen);
end % end: if

% Compute change in fval and individuals in last 'Window' iterations
window = options.stallIterLimit;
if iteration > window
    % The smallest fval in the window should be bestFval.
    % The largest fval in the window should be the oldest one in the
    % window. This value is at iterationIndex+1 (or 1).
    if iterationIdx == window
        % The window runs from index 1:iterationIndex
        maxBestFvalsWindow = bestFvalsWindow(1);
    else
        % The window runs from [iterationIndex+1:end, 1:iterationIndex]
        maxBestFvalsWindow = bestFvalsWindow(iterationIdx+1);
    end % end: if
    funChange = abs(maxBestFvalsWindow-bestFval)/max(1,abs(bestFval));
else
    funChange = Inf;
end % end: if

reasonToStop = '';
exitFlag = [];
if state.iteration >= options.maxGeneration
    reasonToStop = 'Exit: Reaches maximum iteration';
    exitFlag = 0;
elseif toc(state.startTime) > options.maxTime
    reasonToStop = 'Exit: Reaches maximum time';
    exitFlag = -5;
elseif (toc(state.startTime)-state.lastImprovementTime) > options.stallTimeLimit
    reasonToStop = 'Exit: Reaches stall time limit';
    exitFlag = -4;
elseif bestFval < options.objectiveLimit
    reasonToStop = 'Exit: Reaches objective function value limit';
    exitFlag = -3;
elseif funChange <= options.tolFunValue
    reasonToStop = 'Exit: Reaches function value change tolerance';
    exitFlag = 1;
end % end: if

if ~isempty(reasonToStop) && options.verbosity > 0
    fprintf('%s\n', reasonToStop);
    return
end % end: if

% Print header again
if options.verbosity > 1 && ...
    rem(iteration, options.displaySectionSize*options.displayInterval)==0 && ...
    iteration > 0
    fprintf('\n                                 Best            Mean     Stall\n');
    fprintf(  'Iteration     f-count            f(x)            f(x)    Iterations\n');
end % end: if

end % end: function StopParticleswarm



%%%% Function "MergeOptions"

function output = MergeOptions(default, user, name)
%
% Merge a default options struct with a user-defined options struct. Works
% recursively, and will issue warning messages if the user attempts to
% define a field that is not in the default options.
%
% DESCRIPTION:
%
% - All fields in DEFAULT will be present in OUTPUT
% - If a field is in both DEFAULT and USER, then the value from USER is
% present in OUTPUT
% - If a field is present in USER, but not DEFAULT, then issue a warning.
% - Applies recursively
%
% NOTES:
%
%   The argument "name" is optional, and contains a string specifying the
%   name of the options struct. This is primarily used for printing
%   warnings to the user.
%
%   This function works recursively. For example, if there is a struct
%   inside of a struct, then it will recursively apply this merge.
%

% Start by assuming that the OUTPUT is just the DEFAULT
output = default;

% Check if user define option name
if nargin == 2
    structName = '';
else
    structName = [name '.'];
end % end: if

% Merge user-define options with default ones
if ~isempty(user)
    % Check for any overriding fields in the USER-defined struct
    defaultFields = fieldnames(default);
    for i = 1 : length(defaultFields)
        if isfield(user, defaultFields{i})
            C0 = isstruct(default.(defaultFields{i}));
            C1 = isstruct(user.(defaultFields{i}));
            if C0 && C1         % Both are structs
                output.(defaultFields{i}) = MergeOptions(...
                    default.(defaultFields{i}), ...
                    user.(defaultFields{i}), ...
                    [structName defaultFields{i}]);
            elseif ~C0 && ~C1   % Both are fields
                output.(defaultFields{i}) = user.(defaultFields{i});
            elseif C0 && ~C1    %default is struct, user is a field
                disp(['WARNING: ' structName defaultFields{i} ' should be a struct!']);
            elseif ~C0 && C1    %default is struct, user is a field
                disp(['WARNING: ' structName defaultFields{i} ' should not be a struct!']);
            end % end: if
        end % end: if
    end % end: for
    % Check for any fields in USER that are not in DEFAULT
    userFields = fieldnames(user);
    for i = 1 : length(userFields)
        if ~isfield(default, userFields{i})
            disp(['WARNING: unrecognized option: ' structName userFields{i}]);
        end % end: if
    end % end: for
end % end: if

end % end: function MergeOptions

