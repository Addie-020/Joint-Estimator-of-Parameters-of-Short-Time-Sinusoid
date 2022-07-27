function [xBest, fBset, info, dataLog] = ParticleSwarmOptim(objFun, x0, xLb, xUb, options)

% 
% Intelligent optimization algorithm called Particle Swarm Optimization
% Adopted in global search to get a rough estimation of the optimal solution
% 
% Input arguments:
%   @objFun : Object function to be optimized
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

%% Preparation

% Input Vector Size Validation
% ---------------------------
% x0, xLb, xUb: N*1 matrix
% ---------------------------
[n, m] = size(x0);
N = n;
if m ~= 1
    error('x0 is not a column vector!')
end
[n, m] = size(xLb);
if (n ~= N) || (m ~= 1)
    error('xLb is not a valid size!')
end
[n, m] = size(xUb);
if (n ~= N) || (m ~= 1)
    error('xUb is not a valid size!')
end

% Option Defult Set
default.alpha = 0.6;                % Inertia coefficient
default.beta = 0.9;                 % Cognetive coefficient
default.gamma = 0.9;                % Social coefficient
default.nPopulation = 3 * N;        % Population size
default.maxGene = 100;              % Maximum number of generatons
default.tolFun = 1e-6;              % Exit when variance in obejective < tolFun
default.tolX = 1e-10;               % Exit when norm of variance < tolX
default.flagVectorize = false;      % Is the objective function vectorized
default.flagMinimize = true;        % True for minimization, false for maximization
default.xDelMax = xUb - xLb;        % Maximum position update
default.flagWarmStart = false;      % Whether directly use the initial point
default.guessWeight = 0.2;          % On range [0,0.9); 0 for ignore guess, 1 for start at guess
default.plotFun = [];               % Handle a function for plotting the progress
default.display = 'iter';           % Print iteration progress out on the screen
default.printMod = 1;               % Print out every [printMod] iterations

% Set options according to whether user provides options
if nargin == 5
    options = MergeOptions(default, options);
else
    options = default;
end

% Set Optimization Problem: Maximization or Minimization
if options.flagMinimize
    optFun = @min;
else
    optFun = @max;
end

% Check if user provided x0
if isempty(x0)
    x0 = 0.5 * xLb + 0.5 * xUb;
    options.guessWeight = 0.0;
    optFun.flagWarmStart = false;
end


%% Initialization

% ------------------------
% X0, X1, X2: N*M matrix
% X, V: N*M matrix
% ------------------------

% Sample two random points in the search space for each particle
M = options.nPopulation;
X1 = xLb * ones(1, M) + ((xUb - xLb) * ones(1, M)) .* rand(N, M);
X2 = xLb * ones(1, M) + ((xUb - xLb) * ones(1, M)) .* rand(N, M);

% Move initial points towards initial guess, by convex combination
w = options.guessWeight;
X0 = x0 * ones(1, M);
X1 = w * X0 + (1 - w) * X1;
X2 = w * X0 + (1 - w) * X2;

% Initialize population position and velocity
X = X1;         % Initial position of the population
V = X2 - X1;    % Initial velocity of the population

% Check if chosen warm start
if options.flagWarmStart
    X(:, 1) = x0;
    V(:, 1) = zeros(size(x0));
end

% ---------------------------
% F: 1*M matrix
% F_Best: a single number
% P_Best: N*M matrix
% G_Best: N*1 matrix
% ---------------------------

% Calculate Fitness Funtion
if options.flagVectorize
    X_Lb = xLb * ones(1, M);
    X_Ub = xUb * ones(1, M);
    F = objFun(X);
else
    F = zeros(1, M);
    for i = 1 : M
        F(1, i) = objFun(X(:, i));
    end
end

% Find particle best and global best
P_Best = X;                             % Particle best point
F_Best = F;                             % Value of particle best point
[F_Global, G_Idx] = optFun(F_Best);     % Value of best point ever, over all points
G_Best = X(:, G_Idx);                   % Global best point


%% Memory Allocation

% Allocate memory for the dataLog
maxIter = options.maxIter;
dataLog(maxIter) = makeStruct(X, V, F, P_Best, F_Best, G_Best, F_Global, G_Idx);

% Allocate memory for info
info.G_Best = zeros(N, maxIter);
info.F_Global = zeros(1, maxIter);
info.G_Idx = zeros(1, maxIter);
info.P_Best_Var = zeros(N, maxIter);
info.F_Best_Var = zeros(1, maxIter);
info.P_Best_Mean = zeros(N, maxIter);
info.F_Best_Mean = zeros(1, maxIter);
info.P_Var = zeros(N, maxIter);
info.F_Var = zeros(1, maxIter);
info.P_Mean = zeros(N, maxIter);
info.F_Mean = zeros(1, maxIter);
info.iter = 1 : maxIter;


%% Main Loop

% Set parameters
info.exitFlag = 1;
iter = 1;
omega = options.alpha;
c1 = options.beta;
c2 = options.gamma;

% iIteration
while iter <= maxIter
    % Compute new generation of points
    if iter > 1
        r1 = rand(N, M);
        r2 = rand(N, M);
        % Operate according to whether the function is vectorized
        if options.flagVectorize
            V = ...                                     % Update particle velocity
                omega * V + ...                           % Inertia component
                c1 * r1 .* (P_Best - X) + ...             % Cognitive component
                c2 * r2 .* (G_Best * ones(1, M) - X);     % Social component
            X_New = X + V;                              % Update particle postion
            X = max(min(X_New, X_Ub), X_Lb);            % Clamp position to bounds
            F = objFun(X);                              % Calculate fitness function value
            F_Best_New = optFun(F_Best, F);             % Find best fitness value for each
            idxUpdate = (F_Best_New ~= F_Best);         % Index of best points to be updated
            P_Best(:, idxUpdate) = X(:, idxUpdate);     % Update particle best point
            F_Best = F_Best_New;                        % Update value of particle best point
            [F_Global, G_Idx] = optFun(F_Best);         % Update value of best point ever, over all points
            G_Best = X(:, G_Idx);                       % Update global best point
        else
            for idx = 1 : M
                V(:, idx) = ...                                                 % Update particle velocity
                    omega * V(:, idx) + ...                                             % Inertia component
                    c1 * r1(:, idx) .* (P_Best(:, idx) - X(:, idx)) + ...               % Cognitive component
                    c2 * r2(:, idx) .* (G_Best - X(:, idx));                            % Social component
                X_New = X(:, idx) + V(: ,idx);                                  % Update particle position
                X(:, idx) = max(min(X_New, xUpp), xLow);                        % Clamp position to bounds
                F(:, idx) = objFun(X(:, idx));                                  % Calculate fitness function value
                [F_Best(1, idx), iBest] = optFun([F(1, idx), F_Best(1, idx)]);  % Find best fitness value
                if iBest == 1       %Then new point is better!
                    P_Best(:, idx) = X(:, idx);
                    [F_Global, iBest] = optFun([F_Best(1,idx), F_Global]);
                    if iBest == 1   %Then new point is the global best!
                        G_Best = P_Best(:, idx);
                    end
                end
            end
        end
    end

    % Log Data
    dataLog(iter) = makeStruct(X, V, F, P_Best, F_Best, P_Best, F_Global, I_Global);
    info.G_Best(:, iter) = G_Best;
    info.F_Global(iter) = F_Global;
    info.G_Idx(iter) = I_Global;
    info.X_Var(:, iter) = var(X, 0, 2);
    info.P_Best_Var(:, iter) = var(P_Best, 0, 2);
    info.X_Mean(:, iter) = mean(X, 2);
    info.P_Best_Mean(:, iter) = mean(P_Best, 2);
    info.F_Var(1, iter) = var(F);
    info.F_Best_Var(1, iter) = var(F_Best);
    info.F_Mean(1, iter) = mean(F);
    info.F_Best_Mean(1, iter) = mean(F_Best);
    % Plot
    if ~isempty(options.plotFun)
        options.plotFun(dataLog(iter), iter);
    end
    % Print
    xVar = norm(info.X_Var(:, iter));
    if strcmp('iter', options.display)
        if mod(iter - 1, options.printMod)==0
            fprintf('iter: %3d,  fBest: %9.3e,  fVar: %9.3e  xVar: %9.3e  \n',...
                iter, info.F_Global(iter), info.F_Var(1, iter), xVar);
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
end


%% Output
xBest = info.G_Best(:, end);
fBset = info.F_Global(end);
info.input = MakeStruct(objFun, x0, xLb, xUb, options);
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





