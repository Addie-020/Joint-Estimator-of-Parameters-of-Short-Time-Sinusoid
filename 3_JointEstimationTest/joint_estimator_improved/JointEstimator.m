function [xBest, yBest, info] = JointEstimator(xn, Fs, options)

%
% Joint estimator of frequency and phase of sinusoid
% Correlation-based method
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
% Author: Zhiyu Shen @Nanjing University
% Date  : Aug 3, 2022
%

%%% Preparation

% Input Vector Size Validation
n = size(xn, 1);
if n ~= 1
    error('Input signal must be in a row vector!');
end

% Option Defult Set
default.maxIter         = 100;          % Maximum iteration times
default.display         = 3;            % Print iteration progress out on the screen
default.printMod        = 1;            % Print out every [printMod] iterations
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
    optPso.printMod = 'iter';
    optCg.printMod = [];
elseif options.display == 2
    optPso.printMod = [];
    optCg.printMod = 'iter';
else
    optPso.printMod = [];
    optCg.printMod = [];
end


%%% Compute Sequence Information

% Compute mean and variance of test signal
Ns = length(xn);
miu0 = sum(xn) / Ns;
sigma0 = sqrt(sum((xn - repmat(miu0, 1, Ns)).^2) / Ns);

% Compute signal information for correlation computation
Ct = (xn - repmat(miu0, 1, Ns)) ./ repmat(sigma0, 1, Ns);


%%% Memory Allocation

% Allocate memory for info
info.fIter          = zeros(1, maxIter);        % Best frequency value of current iteration
info.pIter          = zeros(1, maxIter);        % Best phase value of current iteration
info.yIter          = zeros(1, maxIter);        % Optimal objective function value of current iteration
info.fGrad          = zeros(1, maxIter);        % Frequency component of gradient of current iteration
info.pGrad          = zeros(1, maxIter);        % Phase component of gradient of current iteration
info.iter           = 1 : maxIter;


%%% Search process

xBest = zeros(1, 2);
yBest = 3;
for iter = 1 : maxIter

    % Global search with random start
    f0 = randi(100) / 100;
    p0 = randi([0 200]) / 100 * pi;
    x0 = [f0; p0];
    xLb = [0; 0];
    xUb = [1; 2*pi];
    [xGlob, ~, ~, ~] = ParticleSwarmOptim(Ct, Fs, x0, xLb, xUb, optPso);

    % Local search
    [xIter, yIter, infoLoc] = ConjGradeOptim(xGlob, Ct, Fs, optCg);

    % Log Data
    info.fIter(iter) = xIter(1);
    info.pIter(iter) = xIter(2);
    info.yIter(iter) = yIter;
    info.fGrad(iter) = infoLoc.gradFreq(end);
    info.pGrad(iter) = infoLoc.gradPha(end);
    
    % Whether new iteration is better
    if yIter < yBest
        xBest = xIter;
        yBest = yIter;
    end

    % Print
    if options.display == 0
        if mod(iter - 1, options.printMod) == 0
            fprintf(['iter: %3d,  freq: %9.3e,  pha: %9.3e  objFun: %9.3e  ' ...
                'freqGrad: %9.3e  phaGrad: %9.3e\n'],...
                iter, info.fIter(iter), info.pIter(iter), info.yIter(iter), ...
                info.fGrad(iter), info.pGrad(iter));
        end
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