% Description:  Test Program for Joint Estimator for MSE Measurement
%               (Single Run)
% Projet:       Joint Estimatior of Frequency and Phase
% Date:         July 28, 2022
% Author:       Zhiyu Shen

clear
close all
clc


%%% Set Up Estimation Options

% Set frequency range
fprintf('Choose an estimation frequency range\n');
while 1
    fprintf(['Use preconfigured range (0,1) or set manually? ' ...
        '(0: preconfigured, 1: manual setting)\n']);
    fConfig = input('Your choise [default 0]:');
    if isempty(fConfig) || (fConfig == 0)
        fLb = 0;
        fUb = 1;
        break;
    elseif fConfig == 1
        fLb = input('Frequency lower bound (Hz):');
        fUb = input('Frequency upper bound (Hz):');
        break;
    else
        fprintf('Invalid input! Type your choise again.\n');
    end
end

% Set phase range
fprintf('Choose an estimation phase range\n');
while 1
    fprintf(['Use preconfigured range (0,2pi) or set manually? ' ...
        '(0: preconfigured, 1: manual setting)\n']);
    pConfig = input('Your choise [default 0]:');
    if isempty(pConfig) || (pConfig == 0)
        pLb = 0;
        pUb = 2*pi;
        break;
    elseif pConfig == 1
        pLb = input('Frequency lower bound (Hz):');
        pUb = input('Frequency upper bound (Hz):');
        break;
    else
        fprintf('Invalid input! Type your choise again.\n');
    end
end
paramRange = [fLb, fUb, pLb, pUb];


%%% Generate Signal to Be Estimated

% Set sampling parameters (Hz)
Fs = input('Sampling frequency(Hz) [10]: ');
if isempty(Fs)
    Fs = 10;
end

% Set sampling time
cycles = input('Number of cycles sampled: [0.5]: ');
if isempty(cycles)
    cycles = 0.5;
end
Tt = cycles / ft;                   % Total time of sampling (s)
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
numEst = 100;                           % Number of estimations
options.maxIter = 5;                    % Search times

% Estimate loop
timeTot = zeros(1, numEst);         % Estimation time for each iteration
fe = zeros(1, numEst);              % Estimated frequency of each iteration
pe = zeros(1, numEst);              % Estimated phase of each iteration
for i = 1 : numEst
    tic
    [xBest, ~, ~] = JointEstimator(xn, Fs, options, optionsPso, optionsGrad);
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
freqErr = fe - ft;
freqMse = sum(freqErr.^2) / numEst;
phaErrVec = abs([pe-pt; pe-pt+2*pi; pe-pt-2*pi]);
phaErr = min(phaErrVec);
phaMse = sum(phaErr.^2) / numEst;

%%% Output

fprintf('\n-------- Input Signal --------\n');
fprintf('Frequency: %.3d Hz\n', ft);
fprintf('Phase: %.3d rad\n', pt);

fprintf('\n-------- Error Analysis --------\n');
fprintf('Frequency Square Error: %.3d\n', freqMse);
fprintf('Phase Square Error: %.3d\n', phaMse);

fprintf('\n-------- Time Used --------\n');
fprintf('Mean Time: %.3f s\n', timeMean);
fprintf('Time Variance: %.3f s\n', timeVar);
fprintf('Total Time: %.3f s\n', timeTot);



