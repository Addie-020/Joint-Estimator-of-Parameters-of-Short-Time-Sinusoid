function [freqMse, phaMse, timeMean, timeVar] = JointEstimatorTest2(xn, ...
    ft, pt, Fs, Tt, numEst, maxIter)

%
% Test function of Joint Estimator
% Run estimator for multiple times
% 
% Input arguments:
%   @xn     : Signal to be estimated
%   @ft     : Frequency of signal to be estimated
%   @pt     : Phase of signal to be estimated
%   @Fs     : Sampling rate (Hz)
%   @Tt     : Total time of sampling (s)
%   @numEst : Estimation times for each test
%   @maxIter: Maximum iteration time for each estimation
%
% Output arguments:
%   @freqMse : MSE of frequency estimated
%   @phaMse  : MSE of phase estimated
%   @timeMean: Mean of time of each estimation
%   @timeVar : Variance of time of each estimation
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Oct 3, 2022
%

%%% Estimation Process

% Define estimator options
options.maxIter = maxIter;
options.display = 3;

% Estimate loop
timeTot = zeros(1, numEst);         % Estimation time for each iteration
fe = zeros(1, numEst);              % Estimated frequency of each iteration
pe = zeros(1, numEst);              % Estimated phase of each iteration
for i = 1 : numEst
    tic
    [xBest, ~, ~] = JointEstimator(xn, Fs, options);
    timeTot(i) = toc;
    % Assign results
    fe(i) = xBest(1);
    pe(i) = xBest(2);
end


%%% Process result

% Process estimation time
timeEst = timeTot + Tt;
timeMean = sum(timeEst) ./ numEst;
timeVar = sum((timeEst-timeMean).^2) / numEst;

% Calculate error
freqMse = sum((fe-ft).^2) / numEst;
phaMse = sum((pe-pt).^2) / numEst;

end



%%%% Function "MergeOptions"

function [xBest, yBest, info] = JointEstimator(xn, Fs, options)

%
% Joint estimator of frequency and phase of sinusoid
% Correlation-based method
% Using self-coded PSO algorithms
% Using MATLAB fmincon as gradient method
% 
% Input arguments:
%   @xn     : Signal to be estimated
%   @Fs     : Sampling rate
%   @options: Optimization options, for more details see 'Option Defult
%             Set' in 'Preparation' part
%
% Output arguments:
%   @xBest  : Optimal point (variable)
%   @fBest  : Optimal value of object function
%   @tTot   : Total time of computation
%   @info   : Information of the optimization process
%   @dataLog: Data log of each iteration
%

%%% Preparation

% Input Vector Size Validation
n = size(xn, 1);
if n ~= 1
    error('Input signal must be in a row vector!');
end

% Option Defult Set
default.maxIter         = 100;          % Maximum iteration times
default.display         = 0;            % Print iteration progress out on the screen
default.printMod        = 3;            % Print out every [printMod] iterations
default.maxRuntime      = 1;            % Maximum run time of each estimation (s)

% Set options according to user inputs
if nargin == 3
    options = MergeOptions(default, options);
else
    options = default;
end

% Assign some paramters
maxIter = options.maxIter;

% Display options
% 0: Display each iteration in joint estimator
% 1: Display each iteration in particle swarm optimization
% 2: Display each iteration in conjugate gradient algorithm
if options.display == 1
    optionPso.Display = 'iter';
    optionGradient.Display = 'off';
elseif options.display == 2
    optionPso.Display = 'none';
    optionGradient.Display = 'iter';
else
    optionPso.Display = 'none';
    optionGradient.Display = 'off';
end

% Options of 'fmincon" function
optionGradient.Algorithm = 'interior-point' ;



%%% Compute Sequence Information

% Compute mean and variance of test signal
Ns = length(xn);
miu0 = sum(xn) / Ns;
sigma0 = sqrt(sum((xn-miu0).^2) / Ns);

% Compute signal information for correlation computation
Ct = (xn-miu0) ./ sigma0;


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

xBest = zeros(1, 2);
yBest = 3;
for iter = 1 : maxIter

    startTime = tic;
    fun = @(X)ObjFun(X, Ct, Fs);
    
    % Global search with random start
    xLb = [0, 0];
    xUb = [1, 2*pi];
    nvars = 2;
    [xGlobal, yGlobal] = particleswarm(fun, nvars, xLb, xUb, optionPso);

    % Local search
    fLbLoc = max(0, xGlobal(1)-0.1);
    fUbLoc = min(1, xGlobal(1)+0.1);
    pLbLoc = max(0, xGlobal(2)-pi/50);
    pUbLoc = min(2*pi, xGlobal(2)+pi/50);
    lb = [fLbLoc, pLbLoc];
    ub = [fUbLoc, pUbLoc];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];
    [xIter, yIter, ~, ~, ~, gIter] = fmincon(fun, xGlobal, A, b, Aeq, beq, lb, ub, nonlcon, optionGradient);
    

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

