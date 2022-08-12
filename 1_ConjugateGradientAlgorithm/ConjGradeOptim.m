function [xBest, yBest, info, dataLog] = ConjGradeOptim(x0, tn, fs, options)
 
%
% Conjugate Gradient Algorithm (with Polak-Ribiere method)
% i.e. PR-CG algorithm
% Adopt golden section method as linear search optimization algorithm
% 
% Input arguments:
%   @x0     : Initial value of variables
%   @tn     : Sequence to be estimated
%   @fs     : Sampling rate of input sequence
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
default.epsilon         = 1e-9;         % Exit falg for 2D search
default.stepErr         = 1e-4;         % Exit flag for 1D search
default.stepDist        = 0.5;          % 1D search distance
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
funVal = ObjFun(xVal, tn, fs);

% Calculate partial differential
xDel = diag(h * ones(D, 1));                % D * D matrix
xInc = repmat(xVal, 1, D) + xDel;           % D * D matrix
funValInc = ObjFun(xInc, tn, fs);           % 1 * D matrix
gVal0 = ((funValInc - funVal * ones(1, D)) / h).';
gVal = gVal0;
dVal = zeros(D, 1);


%%% Memory Allocation

% Allocate memory for the dataLog
dataLog(maxIter) = MakeStruct(xVal, funVal, gVal);

% Allocate memory for info
info.freqVal        = zeros(1, maxIter);        % Frequency value of current iteration
info.phaVal         = zeros(1, maxIter);        % Phase value of current iteration
info.ampVal         = zeros(1, maxIter);        % Amplitude value of current iteration
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

    % Update multi-dimension search direction
    dVal = -gVal + bVal .* dVal;

    % Optimize search step with linear search optimization algorithm
    [xVal, funVal] = StepOptim(xVal, dVal, stepErr, stepDist, tn, fs);
    
    % Calculate partial differential
    xInc = repmat(xVal, 1, D) + xDel;
    funValInc = ObjFun(xInc, tn, fs);
    gVal0 = gVal;
    gVal = ((funValInc - funVal * ones(1, D)) / h).';
    
    % Log Data
    dataLog(iter) = MakeStruct(xVal, funVal, gVal);
    info.freqVal(iter)  = xVal(1);
    info.phaVal(iter)   = xVal(2);
    info.ampVal(iter)   = xVal(3);
    info.funVal(iter)   = funVal;
    info.gradFreq(iter) = abs(gVal(1));
    info.gradPha(iter)  = abs(gVal(2));
    info.gradAmp(iter)  = abs(gVal(3));

    % Print
    if strcmp('iter', options.display)
        if mod(iter - 1, options.printMod) == 0
            fprintf(['iter: %3d,  freq: %9.3e  pha: %9.3e  amp: %9.3e  ' ...
                'objFun: %9.3e  freqErr: %9.3e  phaErr: %9.3e  ampErr: %9.3e\n'],...
                iter, info.freqVal(iter), info.phaVal(iter), info.ampVal(iter), ...
                info.funVal(iter), info.gradFreq(iter), info.gradPha(iter), info.gradAmp(iter));
        end
    end
    
    % Update iteration time
    iter = iter + 1;

end

xBest = xVal;
yBest = funVal;

end



%% Linear Search Optimization

function [xBest, yBest] = StepOptim(x0, d0, epsilon, dist, tn, fs)

%
% Linear search optimization algorithm
% Optimize the step of multi-dimension search
% Adopt golden section method
% 
% Input arguments:
%   @x0     : Initial variable value
%   @g0     : Initial gradient value
%   @funVal0: Initial function value
%   @epsilon: Exit flag for search
%   @dist   : Search distance
%   @tn     : Sequence to be estimated
%   @fs     : Sampling rate of input sequence
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
fc = ObjFun(xc, tn, fs);
xd = x0 + d * d0;
fd = ObjFun(xd, tn, fs);


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
        fc = ObjFun(xc, tn, fs);
    elseif fc < fd
        % Update end points
        a = c;
        % Update end middle points
        c = d;
        d = a + w2 * (b - a);
        % Update function values
        fc = fd;
        xd = x0 + d * d0;
        fd = ObjFun(xd, tn, fs);
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
fa = ObjFun(xa, tn, fs);
xb = x0 + b * d0;
fb = ObjFun(xb, tn, fs);

% Compare function value of end points and determine output value
if fa > fb
    yBest = fa;
    xBest = xa;
else
    yBest = fb;
    xBest = xb;
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