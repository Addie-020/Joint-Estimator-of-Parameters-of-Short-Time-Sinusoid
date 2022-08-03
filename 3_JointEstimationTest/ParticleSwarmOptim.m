function [xBest, fBset, info, dataLog] = ParticleSwarmOptim(Ct, Fs, x0, xLb, xUb, options)

%
% Intelligent optimization algorithm called Particle Swarm Optimization
% Adopted in global search to get a rough estimation of the optimal solution
%
% Input arguments:
%   @Ct     : Necessary information of sequence to be estimated
%   @Fs     : Sampling rate
%   @x0     : Initial value of variables
%   @xLb    : Lower bound of variables
%   @xUb    : Upper bound of variables
%   @options: Optimization options, for more details see 'Option Defult
%             Set' in 'Preparation' part
%
% Output arguments:
%   @xBest  : Optimal point (variable)
%   @fBest  : Optimal value of object function
%   @info   : Information of the optimization process
%   @dataLog:
%
% Author: Zhiyu Shen @Nanjing University
% Date  : July 27, 2022
%

%%% Preparation

% Input Vector Size Validation
% ---------------------------
% x0, xLb, xUb: N*1 matrix
% ---------------------------
[n, m] = size(x0);
D = n;
if m ~= 1
    error('x0 is not a column vector!')
end
[n, m] = size(xLb);
if (n ~= D) || (m ~= 1)
    error('xLb is not of a valid size!')
end
[n, m] = size(xUb);
if (n ~= D) || (m ~= 1)
    error('xUb is not of a valid size!')
end

% Option Defult Set
default.alpha           = 0.6;          % Inertia coefficient
default.beta            = 0.9;          % Cognetive coefficient
default.gamma           = 0.9;          % Social coefficient
default.nPopulation     = 3 * D;        % Population size
default.maxGene         = 1000;         % Maximum number of generatons
default.tolFun          = 1e-9;         % Exit when variance in obejective < tolFun
default.tolX            = 1e-9;         % Exit when norm of variance < tolX
default.xDelMax         = xUb - xLb;    % Maximum position update
default.guessWeight     = 0.2;          % On range [0,0.9); 0 for ignore guess, 1 for start at guess
default.plotFun         = [];           % Handle a function for plotting the progress
default.display         = 'iter';       % Print iteration progress out on the screen
default.printMod        = 1;            % Print out every [printMod] iterations

% Set options according to user inputs
if nargin == 6
    options = MergeOptions(default, options);
else
    options = default;
end

% Check if user provided x0
if isempty(x0)
    x0 = 0.5 * xLb + 0.5 * xUb;
    options.guessWeight = 0.0;
    options.flagWarmStart = false;
end


%%% Initialization

% X0, X1, X2: N*M matrix
% X, V: N*M matrix
% D is searching dimension, i.e. the number of variables
% P is the number of populations

% Sample two random points in the search space for each particle
% Get two set of random particles (D dimensions, P populations)
P = options.nPopulation;
X1 = xLb * ones(1, P) + ((xUb - xLb) * ones(1, P)) .* rand(D, P);
X2 = xLb * ones(1, P) + ((xUb - xLb) * ones(1, P)) .* rand(D, P);

% Move initial points towards initial guess, by convex combination
omega = options.guessWeight;
X0 = x0 * ones(1, P);
X1 = omega * X0 + (1 - omega) * X1;
X2 = omega * X0 + (1 - omega) * X2;

% Initialize population position and velocity
X = X1;         % Initial position of particles
V = X2 - X1;    % Initial velocity of particles

% ---------------------------
% F: 1*P matrix
% F_Best: 1*P matrix
% P_Best: N*M matrix
% G_Best: N*1 matrix
% ---------------------------

% Calculate Fitness Funtion
X_Lb = xLb * ones(1, P);
X_Ub = xUb * ones(1, P);
F = ObjFun(X, Ct, Fs);


% Find particle best and global best
P_Best = X;                             % Particle with best fitness value in each population
F_Best = F;                             % Best fitness value in each population
[F_Glob, G_Idx] = min(F_Best);          % Best fitness value of all particles
G_Best = X(:, G_Idx);                   % Particle with best fitness value of all populations


%%% Memory Allocation

% Allocate memory for the dataLog
maxIter = options.maxGene;
dataLog(maxIter) = MakeStruct(X, V, F, P_Best, F_Best, G_Best, F_Glob, G_Idx);

% Allocate memory for info
info.G_Best         = zeros(D, maxIter);        % Global best particle of current generation
info.F_Glob         = zeros(1, maxIter);        % Fitsness value of global best particle of current generation
info.G_Idx          = zeros(1, maxIter);        % Population index of global best particle of current generation
info.P_Best_Var     = zeros(D, maxIter);        % Variance of best particles from beginning to current generation
info.F_Best_Var     = zeros(1, maxIter);        % Variance of fitness values of best particles from beginning to current generation
info.P_Best_Mean    = zeros(D, maxIter);        % Mean of best particles from beginning to current generation
info.F_Best_Mean    = zeros(1, maxIter);        % Mean of fitness values of best particles from beginning to current generation
info.P_Mean         = zeros(D, maxIter);        % Mean of particles of current generation
info.P_Var          = zeros(D, maxIter);        % Variance of particles of current generation
info.F_Mean         = zeros(1, maxIter);        % Mean of fitness values of current generation
info.F_Var          = zeros(1, maxIter);        % Variance of fitness values of current generation
info.iter           = 1 : maxIter;


%%% Main Loop

% Set parameters
info.exitFlag = 1;
iter = 1;
omega = options.alpha;
c1 = options.beta;
c2 = options.gamma;

% Iteration
while iter <= maxIter
    % Compute new generation of points
    if iter > 1
        r1 = rand(D, P);
        r2 = rand(D, P);
        % Operate according to whether the function is vectorized
        V = ...                                     % Update particle velocity
            omega * V + ...                           % Inertia component
            c1 * r1 .* (P_Best - X) + ...             % Cognitive component
            c2 * r2 .* (G_Best * ones(1, P) - X);     % Social component
        X_New = X + V;                              % Update particle postion
        X = max(min(X_New, X_Ub), X_Lb);            % Clamp position to bounds
        F = ObjFun(X, Ct, Fs);                      % Update fitness value of all particles
        F_Best_New = min(F_Best, F);                % Get best fitness value of each population in a new vector
        idxUpdate = (F_Best_New ~= F_Best);         % Index of particles with best fitness value to be updated
        P_Best(:, idxUpdate) = X(:, idxUpdate);     % Update particle with best fitness value of each population
        F_Best = F_Best_New;                        % Update best fitness value of each population
        [F_Glob, G_Idx] = min(F_Best);              % Update best fitness value of all particles
        G_Best = X(:, G_Idx);                       % Update particle with best fitness value of all populations
    end

    % Log Data
    dataLog(iter) = MakeStruct(X, V, F, P_Best, F_Best, G_Best, F_Glob, G_Idx);
    info.G_Best(:, iter)        = G_Best;
    info.F_Glob(iter)           = F_Glob;
    info.G_Idx(iter)            = G_Idx;
    info.P_Mean(:, iter)        = mean(X, 2);
    info.P_Var(:, iter)         = var(X, 0, 2);
    info.P_Best_Mean(:, iter)   = mean(P_Best, 2);
    info.P_Best_Var(:, iter)    = var(P_Best, 0, 2);
    info.F_Mean(1, iter)        = mean(F);
    info.F_Var(1, iter)         = var(F);
    info.F_Best_Mean(1, iter)   = mean(F_Best);
    info.F_Best_Var(1, iter)    = var(F_Best);

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
if strcmp('iter', options.display) || strcmp('final', options.display)
    switch info.exitFlag
        case 0
            fprintf('PSO Algorithm Converged. Exit: fVar < tolFun\n');
        case 1
            fprintf('Maximum Iteration Reached.\n')
        case 2
            fprintf('Optimization Converged. Exit: norm(xVar) < tolX\n');
    end
end

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