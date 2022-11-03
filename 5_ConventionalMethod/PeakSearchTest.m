% Description:  Test Program for DTFT Peak Search Estimator for Single Run
% Projet:       Joint Estimatior of Frequency and Phase
% Date:         Oct 31, 2022
% Author:       Zhiyu Shen

clear
close all
clc

%%% Generate Signal to Be Estimated

% Set parameter type and define signal to be estimated
fprintf('Set parameter type: fixed, random, user input');
paramType = input('Input type: (f/r/u) [f]:', 's');
if isempty(paramType) || (paramType == 'f')
    ft = 0.01;                              % Frequency of test signal (Hz)
    pt = pi/3;                              % Phase of test signal (rad)
elseif paramType == 'r'
    ft = randi([1 100]) / 100;              % Frequency of test signal (Hz)
    pt = randi([0 200]) * pi / 100;         % Phase of test signal (rad)
elseif paramType == 'u'
    ft = input('Frequency (Hz): ');
    pt = input('Initial phase (rad): ');
else
    error('Invalid input!');
end

% Set sampling parameters (Hz)
Fs = input('Sampling frequency(Hz) [10]: ');
if isempty(Fs)
    Fs = 10;
end
Tt = 5 / ft;                        % Total time of sampling (s)
Ns = round(Tt * Fs);                % Total sampling points

% Generate original signal sequence
xt = (0 : Ns - 1) / Fs;             % Time index
at = 1;                             % Signal amplitude
xn0 = at * cos(2*pi*ft*xt + pt);    % Test signal

% Add noise to signal
addNoise = input('Add noise to signal? Y/N [N]: ', 's');
if isempty(addNoise) || (addNoise == 'N')
    xn = xn0;
elseif addNoise == 'Y'
    % Define SNR
    snrSig = input('SNR(dB) [40]: ');
    if isempty(snrSig)
        snrSig = 40;
    end
    sigmaN = at / 10.^(snrSig/20);      % Standard variance of noise
    sigNoise = sigmaN * randn(1, Ns);   % Additive white Gaussian noise
    xn = xn0 + sigNoise;
end

%%% Estimation Process

% Define estimator options
numEst = 10;                             % Estimation times

% Estimate loop
timeTot = zeros(1, numEst);         % Estimation time for each iteration
fe = zeros(1, numEst);              % Estimated frequency of each iteration
pe = zeros(1, numEst);              % Estimated phase of each iteration
for i = 1 : numEst
    tic
%     [xBest, yBest] = PeakSearchEstimator2(xn, Fs);
    [xBest, yBest] = PhaseDiff(xn, Fs);
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
phaMse = sum((pe-pt).^2) ./ numEst;


fprintf('\n-------- Input Signal --------\n');
fprintf('Frequency: %.3d Hz\n', ft);
fprintf('Phase: %.3d rad\n', pt);

fprintf('\n-------- Error Analysis --------\n');
fprintf('Frequency MSE: %.3d\n', freqMse);
fprintf('Phase MSE: %.3d\n', phaMse);

fprintf('\n-------- Time Used --------\n');
fprintf('Sampling Time: %.3f s\n', Tt);
fprintf('Mean Estimation Time: %.3f s\n', timeMean);



