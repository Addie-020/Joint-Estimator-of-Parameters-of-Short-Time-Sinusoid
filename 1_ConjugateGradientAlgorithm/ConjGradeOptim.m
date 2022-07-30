function [xBest, yBest, info, dataLog] = ConjGradeOptim(x0, an, objFun, options)
 
%
% Conjugate Gradient Algorithm (with Polak-Ribiere method)
% i.e. PR-CG algorithm
% 
% Input arguments:
%   @objFun : Object function to be optimized, must be vectorized
%   @x0     : Initial value of variables
%   @an     : Sequence to be estimated
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
    error('x0 is not a column vector!');
end
[n, m] = size(xLb);
if (n ~= D) || (m ~= 1)
    error('xLb is not a valid size!');
end
[n, m] = size(xUb);
if (n ~= D) || (m ~= 1)
    error('xUb is not a valid size!');
end

% Option Defult Set
default.alpha           = 0.6;          % Step length for 2D search
default.delta           = 1e-9;         % Infinitesimal when calculating gradient
default.epsilon         = 1e-9;         % Exit falg for 2D search
default.stepError       = 1e-9;         % Exit flag for 1D search
default.maxIter         = 100;          % Maximum iteration times
default.display         = 'iter';       % Print iteration progress out on the screen
default.printMod        = 1;            % Print out every [printMod] iterations

% Set options according to user inputs
if nargin == 3
    options = MergeOptions(default, options);
else
    options = default;
end

% Assign some paramters
df = options.delta;
dp = options.delta;
epsilon = options.epsilon;
maxIter = options.maxIter;


%%% Process the initial point

% Calculate function value of initial point
freqVal = x0(1);
phaVal = x0(2);
funVal0 = objFun(x0, an, fs);

% Calculate partial differential on frequency
xIncFreq = [freqVal + df, phaVal];
funValdf = objFun(xIncFreq, an, fs);
gFreq = (funValdf - funVal0) / df;

% Calculate partial differential on phase
xIncPha = [freqVal, phaVal + dp];
funValdp = objFun(xIncPha, an, fs);
gPha = (funValdp - funVal0) / dp;

% Complete differential
g0 = [gFreq, gPha];

end