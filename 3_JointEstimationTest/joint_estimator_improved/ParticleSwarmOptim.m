function [xBest, fBset, info, dataLog] = ParticleSwarmOptim(Ct, Fs, nVars, xLb, xUb, options)

%
% Intelligent optimization algorithm called Particle Swarm Optimization
% Adopted in global search to get a rough estimation of the optimal solution
%
% Input arguments:
%   @Ct     : Necessary information of sequence to be estimated
%   @Fs     : Sampling rate
%   @nVar   : Variable dimension
%   @xLb    : Lower bound of variables
%   @xUb    : Upper bound of variables
%   @options: Optimization options, for more details see 'Option Defult
%             Set' in 'Preparation' part
%
% Output arguments:
%   @xBest  : Optimal point (variable)
%   @fBest  : Optimal value of object function
%   @info   : Information of the optimization process
%   @dataLog: Data log of each iteration
%
% Author: Zhiyu Shen @Nanjing University
% Date  : July 27, 2022
%

%%% Preparation

% Input Vector Size Validation
% ---------------------------
% xLb, xUb: nV*1 matrix
% ---------------------------
nV = nVars;
[n, m] = size(xLb);
if (n ~= nV) || (m ~= 1)
    error('xLb is not of a valid size!')
end
[n, m] = size(xUb);
if (n ~= nV) || (m ~= 1)
    error('xUb is not of a valid size!')
end

% Option Defult Set
default.inertiaRange     = [0.1, 1.1];    % Inertia range
default.selfAdjustment   = 1.49;          % Self adjustment weight
default.socialAdjustment = 1.49;          % Social adjustment weight
default.nParticles       = 3*nV;          % Number of particles
default.minNeighborFrac  = 0.25;          % Minimum neighborhood size fraction
default.maxGene          = 1000;          % Maximum number of generatons
default.initialSwarm     = [];            % Initial particle swarm
default.initialSwarmSpan = 1000;          % Initial swarm span
default.tolFun           = 1e-16;         % Exit when variance in obejective < tolFun
default.tolX             = 1e-16;         % Exit when norm of variance < tolX
default.xDelMax          = xUb - xLb;     % Maximum position update
default.guessWeight      = 0.2;           % On range [0,0.9); 0 for ignore guess, 1 for start at guess
default.plotFun          = [];            % Handle a function for plotting the progress
default.display          = 'iter';        % Print iteration progress out on the screen
default.printMod         = 1;             % Print out every [printMod] iterations

% Set options according to user inputs
if nargin == 6
    options = MergeOptions(default, options);
else
    options = default;
end

% Assign parameters
nP = options.nPopulation;                % Particle(population) size
cSelf = options.selfAdjustment;          % Self adjustment weight
cSocial = options.socialAdjustment;      % Social adjustment weight
minNeighborSize = max(2, ...
    floor(nP*options.minNeighborFrac));  % Minimum neighborhood size
minInertia = options.inertiaRange(1);    % Minimum inertia
maxInertia = options.inertiaRange(2);    % Maximum inertia
xLbMat = repmat(xLb, 1, nP);             % Particle position lower bound matrix
xUbMat = repmat(xUb, 1, nP);             % Particle position upper bound matrix

% X0, X1, X2: N*M matrix
% X, V: N*M matrix
% D is searching dimension, i.e. the number of variables
% P is the number of populations

% Set initial position of particles (D dimensions, nPop populations)
partPos = xLbMat + (xUbMat - xLbMat) .* rand(nV, nP);

% Set initial velocity of particles
partVelo = (xUbMat - xLbMat) .* (rand(nV, nP) * 2 - 1);

% ---------------------------
% F: 1*P matrix
% F_Best: 1*P matrix
% P_Best: N*M matrix
% G_Best: N*1 matrix
% ---------------------------

% Calculate Fitness Funtion
partFval = ObjFun(partPos, Ct, Fs);

% Find particle best and global best
partBest = partPos;                             % Particle with best fitness value in each population
partBestFval = partFval;                        % Best fitness value in each population
[globBestFval, globIdx] = min(partBestFval);    % Best fitness value of all particles
globBest = partPos(:, globIdx);                 % Particle with best fitness value of all populations

% Set iteration coefficients
inertiaCnt = 0;                         % Adaptive inertia counter
inertia = maxInertia;                   % Inertia
nNeighbor = minNeighborhoodSize;        % Adaptive neighborhood size
       

%%% Memory Allocation

% Allocate memory for the dataLog
maxIter = options.maxGene;
dataLog(maxIter) = MakeStruct(partPos, partVelo, partFval, partBest, partBestFval, globBest, globBestFval, globIdx);

% Allocate memory for info
info.globBest         = zeros(nV, maxIter);   % Global best particle of current generation
info.globBestFval     = zeros(1, maxIter);      % Fitsness value of global best particle of current generation
info.globIdx          = zeros(1, maxIter);      % Population index of global best particle of current generation
info.partBestVar      = zeros(nV, maxIter);   % Variance of best particles from beginning to current generation
info.partBestFvalVar  = zeros(1, maxIter);      % Variance of fitness values of best particles from beginning to current generation
info.partBestMean     = zeros(nV, maxIter);   % Mean of best particles from beginning to current generation
info.partBestFvalMean = zeros(1, maxIter);      % Mean of fitness values of best particles from beginning to current generation
info.partMean         = zeros(nV, maxIter);   % Mean of particles of current generation
info.partVar          = zeros(nV, maxIter);   % Variance of particles of current generation
info.partFvalMean     = zeros(1, maxIter);      % Mean of fitness values of current generation
info.partFvalVar      = zeros(1, maxIter);      % Variance of fitness values of current generation
info.iter             = 1 : maxIter;


%%% Iteration

info.exitFlag = 1;
iter = 1;

while iter <= maxIter
    % Generate best neighbor indexs
    bestNeighborIdx = GenBestNeighborIdx(partBestFval, nNeighbor, nP);
    % Compute new generation of points
    if iter > 1
        r1 = rand(nV, nP);
        r2 = rand(nV, nP);
        % Operate according to whether the function is vectorized
        partVelo = ...                                     % Update particle velocity
            omega * partVelo + ...                           % Inertia component
            c1 * r1 .* (partBest - partPos) + ...             % Cognitive component
            c2 * r2 .* (globBest * ones(1, nP) - partPos);     % Social component
        X_New = partPos + partVelo;                              % Update particle postion
        partPos = max(min(X_New, xUbMat), xLbMat);        % Clamp position to bounds
        partFval = ObjFun(partPos, Ct, Fs);                      % Update fitness value of all particles
        F_Best_New = min(partBestFval, partFval);                % Get best fitness value of each population in a new vector
        idxUpdate = (F_Best_New ~= partBestFval);         % Index of particles with best fitness value to be updated
        partBest(:, idxUpdate) = partPos(:, idxUpdate);     % Update particle with best fitness value of each population
        partBestFval = F_Best_New;                        % Update best fitness value of each population
        [globBestFval, globIdx] = min(partBestFval);              % Update best fitness value of all particles
        globBest = partPos(:, globIdx);                       % Update particle with best fitness value of all populations
    end

    % Log Data
    dataLog(iter) = MakeStruct(partPos, partVelo, partFval, partBest, partBestFval, globBest, globBestFval, globIdx);
    info.globBest(:, iter)         = globBest;
    info.globBestFval(iter)        = globBestFval;
    info.globIdx(iter)             = globIdx;
    info.partMean(:, iter)         = mean(partPos, 2);
    info.partVar(:, iter)          = var(partPos, 0, 2);
    info.partBestMean(:, iter)     = mean(partBest, 2);
    info.partBestVar(:, iter)      = var(partBest, 0, 2);
    info.partFvalMean(1, iter)     = mean(partFval);
    info.partFvalVar(1, iter)      = var(partFval);
    info.partBestFvalMean(1, iter) = mean(partBestFval);
    info.partBestFvalVar(1, iter)  = var(partBestFval);

    % Plot
    if ~isempty(options.plotFun)
        options.plotFun(dataLog(iter), iter);
    end

    % Print
    xVar = norm(info.P_Var(:, iter));
    if strcmp('iter', options.display)
        if mod(iter - 1, options.printMod) == 0
            fprintf('iter: %3d,  fBest: %9.3e,  fVar: %9.3e  xVar: %9.3e  \n',...
                iter, info.F_Glob(iter), info.F_Var(1, iter), xVar);
        end
    end

    % Convergence
    if info.F_Var(1, iter) < options.tolFun
        info.exitFlag = 0;
        dataLog = dataLog(1 : iter);
        info = TruncateInfo(info, maxIter, iter);
        break
    elseif xVar < options.tolX
        info.exitFlag = 2;
        dataLog = dataLog(1 : iter);
        info = TruncateInfo(info, maxIter, iter);
        break
    end
    
    % Update iteration time
    iter = iter + 1;

end


%%% Output

xBest = info.G_Best(:, end);
fBset = info.F_Glob(end);
info.input = MakeStruct(Ct, Fs, x0, xLb, xUb, options);
info.fEvalCount = iter * m;

% Print
if strcmp('final', options.display)
    switch info.exitFlag
        case 0
            fprintf('PSO Algorithm Converged. Exit: fVar < tolFun\n');
        case 1
            fprintf('Maximum Iteration Reached.\n')
        case 2
            fprintf('Optimization Converged. Exit: norm(xVar) < tolX\n');
    end
end

end % End of function ParticleSwarmOptim



%%%% Function "MakeState"

function state = MakeState(nvars, xLb, xUb, objFcn, options)
% 
% Create an initial set of particles and objective function values
% 

% A variety of data used in various places
state = struct;
state.iteration = 0;                % Current generation counter
state.startTime = tic;              % Tic identifier
state.stopFlag = false;             % OutputFcns flag to end the optimization
state.lastImprovement = 1;          % Generation stall counter
state.sastImprovementTime = 0;      % Stall time counter
state.funEval = 0;                  % Function value
numParticle = options.SwarmSize;    % Number of particles

% If InitialSwarm is partly empty use the creation function to generate
% population (CreationFcn can utilize InitialSwarm)
if numParticles ~= size(options.InitialSwarm,1)
    state.Positions = feval(options.CreationFcn,problemStruct);
else % the initial swarm was passed in
    state.Positions = options.initialSwarm;
end

% Enforce bounds
if any(any(state.Positions < lbMatrix)) || any(any(state.Positions > ubMatrix))
    state.Positions = max(lbMatrix, state.Positions);
    state.Positions = min(ubMatrix, state.Positions);
end

% Initialize velocities by randomly sampling over the smaller of initSwarmSpan or ub-lb
% min will be initSwarmSpan if either lb or ub is not finite
vMax = min(xUb-xLb, options.initialSwarmSpan);
state.velocities = repmat(-vMax, numParticle, 1) + ...
    repmat(2*vMax, numParticle, 1) .* rand(numParticle, nvars);

% Calculate the objective function for all particles.
% Vectorized call to objFcn
fvals = objFcn(state.Positions);
state.Fvals = fvals(:);
state.FunEval = numParticle;

state.IndividualBestFvals = state.Fvals;
state.IndividualBestPositions = state.Positions;
end % End of function MakeState


%%%% Function "GenBestNeighborIdx"

function bestNeighborIdx = GenBestNeighborIdx(partBestFval, nNeighbor, nPop)
% 
% Generate best neighborhood index
% The best particle in random neighborhood
% The size is controlled by the adaptiveNeighborhoodSize parameter
% 

neighborIdx = zeros(nPop, nNeighbor);
neighborIdx(:, 1) = 1 : nPop;              % First neighbor is self
for i = 1 : nPop
    % Determine random neighbors that exclude the particle itself,
    % which is (nNeighbor-1) particles
    neighbors = randperm(nPop-1, nNeighbor-1);
    % Add 1 to indicies that are >= current particle index
    iShift = neighbors >= i;
    neighbors(iShift) = neighbors(iShift) + 1;
    neighborIdx(i, 2:end) = neighbors;
end

% Identify the best neighbor
[~, bestRowIdx] = min(partBestFval(neighborIdx), [], 2);

% Create the linear index into neighborIndex
bestLinearIdx = (bestRowIdx.'-1).*nPop + (1:nPop);
bestNeighborIdx = neighborIdx(bestLinearIdx);

end



%%%% Function "genBestNeighborIdx"

function newVelo = updateVelo(partVelo, inertia, ...
    bestNeighborIdx, cSelf, cSocial, pIdx, nVar)
% 
% Update the velocities of particles with indices pIdx
% 

% Generate random number distributions for self and social components
randSelf = rand(numel(pIdx), nVar);
randSocial = rand(numel(pIdx), nVar);

oldVelo = partVelo(pIdx,:);

% Update rule
newVelocities = inertia*oldVelo + ...
    cSelf*randSelf.*(state.IndividualBestPositions(pIdx,:)-state.Positions(pIdx,:)) + ...
    cSocial*randSocial.*(state.IndividualBestPositions(bestNeighborIdx(pIdx), :)-state.Positions(pIdx,:));

% Ignore infinite velocities
tfInvalid = ~all(isfinite(newVelocities), 2);
newVelo(tfInvalid) = oldVelo(tfInvalid);

end



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
end

% Merge user-define options with default ones
if ~isempty(user)
    % Check for any overriding fields in the USER-defined struct
    default_fields = fieldnames(default);
    for i = 1 : length(default_fields)
        if isfield(user, default_fields{i})
            C0 = isstruct(default.(default_fields{i}));
            C1 = isstruct(user.(default_fields{i}));
            if C0 && C1         % Both are structs
                output.(default_fields{i}) = MergeOptions(...
                    default.(default_fields{i}), ...
                    user.(default_fields{i}), ...
                    [structName default_fields{i}]);
            elseif ~C0 && ~C1   % Both are fields
                output.(default_fields{i}) = user.(default_fields{i});
            elseif C0 && ~C1    %default is struct, user is a field
                disp(['WARNING: ' structName default_fields{i} ' should be a struct!']);
            elseif ~C0 && C1    %default is struct, user is a field
                disp(['WARNING: ' structName default_fields{i} ' should not be a struct!']);
            end
        end
    end

    % Check for any fields in USER that are not in DEFAULT
    user_fields = fieldnames(user);
    for i = 1 : length(user_fields)
        if ~isfield(default, user_fields{i})
            disp(['WARNING: unrecognized option: ' structName user_fields{i}]);
        end
    end

end

end



%%%% Function "MakeStruct"

function S = MakeStruct(varargin)
%
% A struct is created with the property that each field corresponds to one
% of the arguments passed to this function.
%
% Example:
%
%   If defines:
%       a = 1;
%       b = 2;
%       c = 0;
%       S = makeStruct(a,b,c);
%   Then
%       S.a = 1;
%       S.b = 2;
%       S.c = 0;
%
% Notes:
%
%   Input names should be unique.
%

N_Inputs = length(varargin);

for i = 1 : N_Inputs
    name = inputname(i);
    S.(name) = varargin{i};
end

end



%%%% Function "TruncateInfo"

function info = TruncateInfo(info,maxIter,iter)
%
% Removes the empty entries in the info struct
%

names = fieldnames(info);
for i = 1 : length(names)
    if (isnumeric(info.(names{i})))   % Check if it's a matrix
        if size(info.(names{i}), 2) == maxIter    % Check if it is iteration data
            info.(names{i}) = info.(names{i})(:, 1 : iter);
        end
    end
end

end