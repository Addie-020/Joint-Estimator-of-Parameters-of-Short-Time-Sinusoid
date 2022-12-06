% Description:  Test Program for Joint Estimator for Varying Sampling Time
% Projet:       Joint Estimatior of Frequency and Phase
% Date:         Oct 3, 2022
% Author:       Zhiyu Shen

clear
close all
clc

%% Generate Signal to Be Estimated

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

% Add noise to signal
addNoise = input('Add noise to signal? Y/N [N]: ', 's');
if isempty(addNoise) || (addNoise == 'N')
    noiseFlag = 0;
elseif addNoise == 'Y'
    % Define SNR
    snrSig = input('SNR(dB) [40]: ');
    if isempty(snrSig)
        snrSig = 40;
    end
    noiseFlag = 1;
end


%% Iteration

numCycle = 0.2 : 0.1 : 2.1;         % Number of cycles
Tt = numCycle / ft;                 % Total time of sampling (s)
numTt = length(numCycle);           % Iteration times
freqMse = zeros(1, numTt);          % MSE of frequency
phaMse = zeros(1, numTt);           % MSE of phase
timeMean = zeros(1, numTt);         % Mean of time
timeVar = zeros(1, numTt);          % Variance of time

poolobj = parpool(12);
parfor i = 1 : numTt
    
    Ns = round(Tt(i) * Fs);             % Total sampling points
    
    % Generate original signal sequence
    xt = (0 : Ns - 1) / Fs;             % Time index
    at = 1;                             % Signal amplitude
    xn0 = at * cos(2*pi*ft*xt + pt);    % Test signal

    % Define estimator options
    maxIter = 10;                       % Maximum iteration time for each estimation
    numEst = 50;                        % Estimation times for each test

    % Add noise with varying SNR and estimate
    if ~noiseFlag
        xn = xn0;
    else
        sigmaN = at / 10.^(snrSig/20);      % Standard variance of noise
        sigNoise = sigmaN * randn(1, Ns);   % Additive white Gaussian noise
        xn = xn0 + sigNoise;
    end

    [freqMse(i), phaMse(i), timeMean(i), timeVar(i)] = JointEstimatorTest(xn, ft, pt, Fs, ...
        Tt(i), numEst, maxIter);

    fprintf('Estimation No.%d, Number of cycles = %.1f\n', i, numCycle(i));

end
delete(poolobj);


%% Plot

% Plot relationship between mean estimation time and SNR
timePlt = figure(1);
timePlt.Name = "Relationship between Time and Estimation";
timePlt.WindowState = 'maximized';
% Plot curve
hold on
plot(numCycle, timeMean, 'LineWidth', 2, 'Color', '#D95319', 'Marker', '*', 'MarkerSize', 8);
hold off
% Set the plotting properties
xlabel("Number of Cycles", "Interpreter", "latex");
ylabel("Mean Estimation Time (s)", "Interpreter", "latex");
set(gca, 'Fontsize', 20);

% Plot relationship between MSE and SNR
errPlt = figure(2);
errPlt.Name = "Relationship between MSE and SNR";
errPlt.WindowState = 'maximized';
% Plot frequency MSE-SNR curve
subplot(2, 1, 1);
hold on
plot(numCycle, log10(freqMse), 'LineWidth', 2, 'Color', '#0072BD', 'Marker', '*', 'MarkerSize', 8);
% plot(snrSig, log10(varLb), 'LineWidth', 2, 'Color', '#D95319', 'Marker', 'o', 'MarkerSize', 8);
hold off
xlabel("Number of Cycles", "Interpreter", "latex");
ylabel("$\log_{10}(MSE_{frequency})$", "Interpreter", "latex");
% legend('Joint Estimator', 'CRLB');
set(gca, 'Fontsize', 20);
% Plot phase MSE-SNR curve
subplot(2, 1, 2);
hold on
plot(numCycle, log10(phaMse), 'LineWidth', 2, 'Color', '#D95319', 'Marker', '*', 'MarkerSize', 8);
hold off
xlabel("Number of Cycles", "Interpreter", "latex");
ylabel("$\log_{10}(MSE_{phase})$", "Interpreter", "latex");
set(gca, 'Fontsize', 20);


%% Print Estimation Information

fprintf('\n');
fprintf('Signal frequency: %.3f Hz\n', ft);
fprintf('Signal phase: %.3f rad\n', pt);
if ~noiseFlag
    fprintf('Noise not added.\n');
else
    fprintf('SNR: %.3f dB\n', snrSig);
end



