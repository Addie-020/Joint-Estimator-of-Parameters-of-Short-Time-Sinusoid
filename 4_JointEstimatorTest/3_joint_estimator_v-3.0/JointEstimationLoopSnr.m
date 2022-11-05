% Description:  Test Program for Joint Estimator for Varying SNR
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
Tt = 0.5 / ft;                      % Total time of sampling (s)
Ns = round(Tt*Fs);                  % Total sampling points

% Generate original signal sequence
xt = (0 : Ns-1) / Fs;               % Time index
at = 1;                             % Signal amplitude
xn0 = at * cos(2*pi*ft*xt + pt);    % Test signal

% Define estimator options
maxIter = 10;                       % Maximum iteration time for each estimation
numEst = 1000;                      % Estimation times for each test


%% Add noise with varying SNR and estimate

snrSig = 0 : 5 : 80;                % SNR (dB)
sigmaN = at ./ 10.^(snrSig/20);     % Standard variance of noise
numSnr = length(snrSig);            % Number of different SNRs
freqMse = zeros(1, numSnr);         % MSE of frequency
phaMse = zeros(1, numSnr);          % MSE of phase
timeMean = zeros(1, numSnr);        % Mean of time
timeVar = zeros(1, numSnr);         % Variance of time

poolobj = parpool(12);
parfor i = 1 : numSnr
    
    sigNoise = sigmaN(i)*randn(1, Ns);      % Additive white Gaussian noise
    xn = xn0 + sigNoise;

    [freqMse(i), phaMse(i), timeMean(i), timeVar(i)] = JointEstimatorTest(xn, ft, pt, Fs, ...
        Tt, numEst, maxIter);

    fprintf('Estimation No.%d, SNR = %.1f\n', i, snrSig(i));

end
delete(poolobj);


%% Calculate CRLB
[~, varLbFreq, varLbPha] = CramerRaoCompute(Fs, at, sigmaN, Ns);


%% Plot

% Plot relationship between mean estimation time and SNR
timePlt = figure(1);
timePlt.Name = "Relationship between Time and Estimation";
timePlt.WindowState = 'maximized';
% Plot curve
hold on
plot(snrSig, timeMean, 'LineWidth', 2, 'Color', '#D95319', 'Marker', '*', 'MarkerSize', 8);
hold off
% Set the plotting properties
xlabel("SNR (dB)", "Interpreter", "latex");
ylabel("Mean Estimation Time (s)", "Interpreter", "latex");
ylim([Tt-2, Tt+2]);
set(gca, 'Fontsize', 20);

% Plot relationship between MSE and SNR
errPlt = figure(2);
errPlt.Name = "Relationship between MSE and SNR";
errPlt.WindowState = 'maximized';
% Plot frequency MSE-SNR curve
subplot(2, 1, 1);
hold on
plot(snrSig, log10(freqMse), 'LineWidth', 2, 'Color', '#0072BD', 'Marker', '*', 'MarkerSize', 8);
plot(snrSig, log10(varLbFreq), 'LineWidth', 2, 'Color', '#D95319', 'Marker', 'o', 'MarkerSize', 8);
hold off
xlabel("SNR (dB)", "Interpreter", "latex");
ylabel("$\log_{10}(MSE_{frequency})$", "Interpreter", "latex");
legend('Joint Estimator', 'CRLB');
set(gca, 'Fontsize', 20);
% Plot phase MSE-SNR curve
subplot(2, 1, 2);
hold on
plot(snrSig, log10(phaMse), 'LineWidth', 2, 'Color', '#0072BD', 'Marker', '*', 'MarkerSize', 8);
plot(snrSig, log10(varLbPha), 'LineWidth', 2, 'Color', '#D95319', 'Marker', 'o', 'MarkerSize', 8);
hold off
xlabel("SNR (dB)", "Interpreter", "latex");
ylabel("$\log_{10}(MSE_{phase})$", "Interpreter", "latex");
legend('Joint Estimator', 'CRLB');
set(gca, 'Fontsize', 20);


%% Print Estimation Information

fprintf('\n');
fprintf('Signal frequency: %.3f Hz\n', ft);
fprintf('Signal phase: %.3f rad\n', pt);


