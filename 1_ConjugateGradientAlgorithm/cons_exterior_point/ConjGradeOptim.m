function [xBest, yBest, info] = ConjGradeOptim(x0, Ct, Fs, options)
 
%
% Conjugate Gradient Algorithm (with Polak-Ribiere method)
% i.e. PR-CG algorithm
% Adopt golden section method as linear search optimization algorithm
% 
% Input arguments:
%   @x0     : Initial value of variables
%   @Ct     : Necessary information of sequence to be estimated
%   @Fs     : Sampling rate
%   @options: Optimization options, for more details see 'Option Defult
%             Set' in 'Preparation' part
%
% Output arguments:
%   @xBest  : Optimal point (variable)
%   @fBest  : Optimal value of object function
%   @info   : Information of the optimization process
%
% Author: Zhiyu Shen @Nanjing University
% Date  : July 27, 2022
%

%%% Preparation

% Input Vector Size Validation
% ---------------------------
% x0: D*1 matrix
% ---------------------------
[n, m] = size(x0);
D = n;
if m ~= 1
    error('x0 is not a column vector!');
end

% Option Defult Set
default.alpha           = 0.6;          % Step length for 2D search
default.delta           = 1e-9;         % Infinitesimal for calculating gradient
default.epsilon         = 1e-12;        % Exit falg for 2D search
default.stepErr         = 1e-4;         % Exit flag for 1D search
default.stepDist        = 0.05;         % 1D search distance
default.maxIter         = 100;          % Maximum iteration times
default.display         = 'iter';       % Print iteration progress out on the screen
default.printMod        = 1;            % Print out every [printMod] iterations

% Set options according to user inputs
if nargin == 4
    options = MergeOptions(default, options);
else
    options = default;
end

% Assign some paramters
h = options.delta;
epsilon = options.epsilon;
stepErr = options.stepErr;
stepDist = options.stepDist;
maxIter = options.maxIter;


%%% Process the initial point

% Calculate function value of initial point
xVal = x0;
funVal = ObjFun(xVal, Ct, Fs);

% Calculate partial differential
xDel = diag(h * ones(D, 1));                % D * D matrix
xInc = repmat(xVal, 1, D) + xDel;           % D * D matrix
funValInc = ObjFun(xInc, Ct, Fs);           % 1 * D matrix
gVal0 = ((funValInc - funVal * ones(1, D)) / h).';
gVal = gVal0;
dVal = zeros(D, 1);


%%% Memory Allocation

% Allocate memory for info
info.freqVal        = zeros(1, maxIter);        % Frequency value of current iteration
info.phaVal         = zeros(1, maxIter);        % Phase value of current iteration
info.funVal         = zeros(1, maxIter);        % Objective function value of current iteration
info.gradFreq       = zeros(1, maxIter);        % Frequency component of gradient of current iteration
info.gradPha        = zeros(1, maxIter);        % Phase component of gradient of current iteration
info.iter           = 1 : maxIter;


%%% Main Loop

iter = 1;
while (max(abs(gVal)) > epsilon) && (iter <= maxIter)
    
    % Using PR method to calculate step length of direction updating
    bValTemp = (gVal' * (gVal - gVal0)) / (norm(gVal0))^2;
    bVal = max(bValTemp, 0);

    % Update 2D search direction
    dVal = -gVal + bVal .* dVal;

    % Optimize search step with linear search optimization algorithm
    [xVal, funVal] = StepOptim(xVal, dVal, stepErr, stepDist, Ct, Fs);
    
    % Calculate partial differential
    xInc = repmat(xVal, 1, D) + xDel;
    funValInc = ObjFun(xInc, Ct, Fs);
    gVal0 = gVal;
    gVal = ((funValInc - funVal * ones(1, D)) / h).';
    
    % Log Data
    info.freqVal(iter)  = xVal(1);
    info.phaVal(iter)   = xVal(2);
    info.funVal(iter)   = funVal;
    info.gradFreq(iter) = abs(gVal(1));
    info.gradPha(iter)  = abs(gVal(2));

    % Print
    if strcmp('iter', options.display)
        if mod(iter - 1, options.printMod) == 0
            fprintf(['iter: %3d,  freq: %9.3e,  pha: %9.3e  objFun: %9.3e  ' ...
                'freqGrad: %9.3e  phaGrad: %9.3e\n'],...
                iter, info.freqVal(iter), info.phaVal(iter), info.funVal(iter), ...
                info.gradFreq(iter), info.gradPha(iter));
        end
    end
    
    % Update iteration time
    iter = iter + 1;

end

xBest = xVal;
yBest = funVal;

end



%%%% Linear Search Optimization

function [xBest, yBest] = StepOptim(x0, d0, epsilon, dist, Ct, Fs)

%
% Linear search optimization algorithm
% Optimize the step of 2D search
% Adopt golden section method
% 
% Input arguments:
%   @x0     : Initial variable value
%   @g0     : Initial gradient value
%   @funVal0: Initial function value
%   @epsilon: Exit flag for search
%   @dist   : Search distance
%   @Ct     : Necessary information of sequence to be estimated
%   @Fs     : Sampling rate
%
% Output arguments:
%   @xBest  : Optimal point (variable)
%   @fBest  : Optimal value of object function
%

%%% Initialization

% Calculate golden section
w1 = (sqrt(5) - 1) / 2;     % 0.618
w2 = (3 - sqrt(5)) / 2;     % 0.382

% Set initial scale
a = 0;                      % Upper bound
b = dist;                   % Lower bound

% Set initial middle points
c = a + w1 * (b - a);
d = a + w2 * (b - a);

% Calculate function value of middle points
xc = x0 + c * d0;
fc = ObjFun(xc, Ct, Fs);
xd = x0 + d * d0;
fd = ObjFun(xd, Ct, Fs);


%%% Iteration

iter = 0;
while abs(b - a) > epsilon
    
    % Update points
    if fc > fd
        % Update end points
        b = d;
        % Update end middle points
        d = c;
        c = a + w1 * (b - a);
        % Update function values
        fd = fc;
        xc = x0 + c * d0;
        fc = ObjFun(xc, Ct, Fs);
    elseif fc < fd
        % Update end points
        a = c;
        % Update end middle points
        c = d;
        d = a + w2 * (b - a);
        % Update function values
        fc = fd;
        xd = x0 + d * d0;
        fd = ObjFun(xd, Ct, Fs);
    else
        % Update end points
        a = c;
        b = d;
        % Update end middle points
        c = a + w1 * (b - a);
        d = a + w2 * (b - a);
        % Update function values
        xc = x0 + c * d0;
        fc = ObjFun(xc, Ct, Fs);
        xd = x0 + d * d0;
        fd = ObjFun(xd, Ct, Fs);
    end 

    % Update iteration time
    iter = iter + 1;

end

%%% Set output value

% Calculate function value of endpoints
xa = x0 + a * d0;
fa = ObjFun(xa, Ct, Fs);
xb = x0 + b * d0;
fb = ObjFun(xb, Ct, Fs);

% Compare function value of end points and determine output value
if fa > fb
    yBest = fa;
    xBest = xa;
else
    yBest = fb;
    xBest = xb;
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

