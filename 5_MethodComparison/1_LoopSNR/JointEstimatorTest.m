function [freqMse, phaMse, timeMean, timeVar] = JointEstimatorTest(xn, ...
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
freqMse = sum((fe-ft).^2) ./ numEst;
phaMse = sum((pe-pt-pi/2).^2) ./ numEst;

end
