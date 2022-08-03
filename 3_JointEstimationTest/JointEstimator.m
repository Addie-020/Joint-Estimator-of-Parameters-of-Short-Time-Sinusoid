function [xBest, yBest, info, dataLog] = JointEstimator(xn, Fs, options)

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
%   @info   : Information of the optimization process
%   @dataLog: Data log of each iteration
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Aug 3, 2022
%

%%% Preparation

% Input Vector Size Validation
[n, m] = size(xn);
N = m;
if n ~= 1
    error('Input signal must be in a row vector!');
end

% Option Defult Set
default.maxIter         = 100;          % Maximum iteration times
default.display         = 0;            % Print iteration progress out on the screen
default.printMod        = 1;            % Print out every [printMod] iterations

% Set options according to user inputs
if nargin == 3
    options = MergeOptions(default, options);
else
    options = default;
end

% Assign some paramters
maxIter = options.maxIter;
if options.display == 0
    optPso.printMod = [];
    optCg.printMod = [];
elseif options.display == 1
    optPso.printMod = 'iter';
    optCg.printMod = [];
elseif options.display == 2
    optPso.printMod = [];
    optCg.printMod = 'iter';
end

%%% Compute sequence information

% Compute mean and variance of test signal
miu0 = sum(xn) / Ns;
sigma0 = sqrt(sum((xn - repmat(miu0, 1, Ns)).^2) / Ns);

% Compute signal information for correlation computation
Ct = (xn - repmat(miu0, 1, Ns)) ./ repmat(sigma0, 1, Ns);