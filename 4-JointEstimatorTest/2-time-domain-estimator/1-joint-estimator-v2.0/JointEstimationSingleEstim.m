% Description:  Test Program for Joint Estimator for Single Run
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

% Set parameter type and define signal to be estimated
fprintf('Set parameter type: fixed, random, user input');
paramType = input('Input type: (f/r/u) [f]:', 's');
if isempty(paramType) || (paramType == 'f')
    ft = (fLb+fUb)/2;
    pt = (pLb+pUb)/2;
elseif paramType == 'r'
    ft = fLb + 0.01*randi([0 round(100*(fUb-fLb))]);
    pt = pLb + 0.01*randi([0 round(100*(pUb-pLb))]);
elseif paramType == 'u'
    fInText = ['Frequency (', num2str(fLb), '~', num2str(fUb), ' Hz): '];
    pInText = ['Phase (', num2str(pLb), '~', num2str(pUb), ' rad): '];
    ft = input(fInText);
    pt = input(pInText);
else
    error('Invalid input!');
end

% Set sampling parameters (Hz)
Fs = input('Sampling frequency(Hz) [10*fUb]: ');
if isempty(Fs)
    Fs = fUb*10;
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
M = 10;                             % Search times
options.maxIter = M;

% Estimate loop
tic
[xBest, yBest, info] = JointEstimator(xn, Fs, paramRange, options, [], []);
totTime = toc;

% Assign results
fe = xBest(1);
pe = xBest(2);

% Calculate error
freqErr = (fe-ft).^2;
phiErr = min(abs([pe-pt; pe-pt+2*pi; pe-pt-2*pi])).^2;

fprintf('\n-------- Input Signal --------\n');
fprintf('Frequency: %.3d Hz\n', ft);
fprintf('Phase: %.3d rad\n', pt);

fprintf('\n-------- Estimation Result --------\n');
fprintf('Frequency: %.3d Hz\n', fe);
fprintf('Phase: %.3d rad\n', pe);

fprintf('\n-------- Error Analysis --------\n');
fprintf('Frequency Square Error: %.3d\n', freqErr);
fprintf('Phase Square Error: %.3d\n', phiErr);

fprintf('\n-------- Time Used --------\n');
fprintf('Sampling time: %.3f s\n', Tt);
fprintf('Total time: %.3f s\n', totTime);
fprintf('Mean time: %.3f s\n', info.meanTime);
fprintf('Test time: %.3f s\n', totTime+Tt);



