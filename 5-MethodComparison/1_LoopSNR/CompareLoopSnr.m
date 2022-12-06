% Description:  Comparison of Estimators with Varying SNR
% Projet:       Joint Estimatior of Frequency and Phase
% Date:         Nov 1, 2022
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

% Set sampling cycles
cycles = input('Number of cycles sampled: [0.5]: ');
if isempty(cycles)
    cycles = 0.5;
end

Tt = cycles/ft;                        % Total time of sampling (s)
Tt1 = cycles/ft;
Ns = round(Tt*Fs);                  % Total sampling points
Ns1 = round(Tt1*Fs);

% Generate original signal sequence
xt = (0:Ns-1)/Fs;                   % Time index
xt1 = (0:Ns-1)/Fs;
at = 1;                             % Signal amplitude
xn0 = at*sin(2*pi*ft*xt+pt);        % Test signal
xn1 = at*cos(2*pi*ft*xt1+pt);

% Define estimator options
maxIter = 5;                        % Maximum iteration time for each estimation
numEst = 200;                       % Estimation times for each test


%% Add noise with varying SNR and estimate

snrSig = 20 : 5 : 80;               % SNR (dB)
sigmaN = at ./ 10.^(snrSig/20);     % Standard variance of noise
numSnr = length(snrSig);            % Number of different SNRs
freqMseA = zeros(1, numSnr);        % MSE of frequency (Joint)
phaMseA = zeros(1, numSnr);         % MSE of phase (Joint)
timeMeanA = zeros(1, numSnr);       % Mean of time (Joint)
timeVarA = zeros(1, numSnr);        % Variance of time (Joint)
freqMseB = zeros(1, numSnr);        % MSE of frequency (Improved)
phaMseB = zeros(1, numSnr);         % MSE of phase (Improved)
timeMeanB = zeros(1, numSnr);       % Mean of time (Improved)
timeVarB = zeros(1, numSnr);        % Variance of time (Improved)
freqMseC = zeros(1, numSnr);        % MSE of frequency (Peak)
phaMseC = zeros(1, numSnr);         % MSE of phase (Peak)
timeMeanC = zeros(1, numSnr);       % Mean of time (Peak)
timeVarC = zeros(1, numSnr);        % Variance of time (Peak)
freqMseD = zeros(1, numSnr);        % MSE of frequency (Phase)
phaMseD = zeros(1, numSnr);         % MSE of phase (Phase)
timeMeanD = zeros(1, numSnr);       % Mean of time (Phase)
timeVarD = zeros(1, numSnr);        % Variance of time (Phase)

poolobj = parpool(12);
parfor i = 1 : numSnr
    
    % Add noise
    sigNoise = sigmaN(i)*randn(1, Ns);      % Additive white Gaussian noise
    xn = xn0 + sigNoise;
    xn2 = xn1 + sigNoise;
    
    % Estimate with Joint Estimator
    [freqMseA(i), phaMseA(i), timeMeanA(i), timeVarA(i)] = JointEstimatorTest(xn, ...
        ft, pt, Fs, Tt, numEst, maxIter);
    % Estimate with Improved Joint Estimator
    [freqMseB(i), phaMseB(i), timeMeanB(i), timeVarB(i)] = JointEstimatorTest2(xn, ...
        ft, pt, Fs, Tt, numEst, maxIter)
    % Estimate with DTFT Peak Search
    [freqMseC(i), phaMseC(i), timeMeanC(i), timeVarC(i)] = PeakSearchTest(xn2, ...
        ft, pt, Fs, Tt1, numEst)
    % Estimate with Phase Difference Method
    [freqMseD(i), phaMseD(i), timeMeanD(i), timeVarD(i)] = PhaseDiffTest(xn2, ...
        ft, pt, Fs, Tt1, numEst)


    fprintf('Estimation No.%d, SNR = %.1f\n', i, snrSig(i));

end
delete(poolobj);


%% Calculate CRLB
[~, varLbFreq, varLbPha] = CramerRaoCompute(Fs, at, sigmaN, Ns);


%% Plot

% Text displayed
textParam = ['$F_s=', num2str(Fs), '$ Hz, $t_s=', num2str(cycles), ...
    'T$, $f_{test}=', num2str(ft), '$ Hz, $\phi_{test}=', ...
    num2str(pt/pi), '\pi$ rad'];

% Plot relationship between mean estimation time and SNR
timePlt = figure(1);
timePlt.Name = "Relationship between Time and Estimation";
timePlt.WindowState = 'maximized';
% Plot curve
hold on
plot(snrSig, timeMeanA, 'LineWidth', 2, 'Color', '#4DBEEE', 'Marker', '*');
plot(snrSig, timeMeanB, 'LineWidth', 2, 'Color', '#77AC30', 'Marker', 'o');
plot(snrSig, timeMeanC, 'LineWidth', 2, 'Color', '#D95319', 'Marker', '+');
plot(snrSig, timeMeanD, 'LineWidth', 2, 'Color', '#EDB120', 'Marker', '.');
hold off
% Set the plotting properties
title(textParam, 'Interpreter', 'latex');
xlabel("SNR (dB)", "Interpreter", "latex");
ylabel("Mean Estimation Time (s)", "Interpreter", "latex");
legend('Joint Estimator', 'Improved Joint', 'Peak Search', 'Phase Difference');
set(gca, 'Fontsize', 20);

% Plot relationship between MSE and SNR
errPlt = figure(2);
errPlt.Name = "Relationship between MSE and SNR";
errPlt.WindowState = 'maximized';
% Plot frequency MSE-SNR curve
subplot(2, 1, 1);
hold on
plot(snrSig, log10(freqMseA), 'LineWidth', 2, 'Color', '#4DBEEE', 'Marker', '*');
plot(snrSig, log10(freqMseB), 'LineWidth', 2, 'Color', '#77AC30', 'Marker', 'o');
plot(snrSig, log10(freqMseC), 'LineWidth', 2, 'Color', '#D95319', 'Marker', '+');
plot(snrSig, log10(freqMseD), 'LineWidth', 2, 'Color', '#EDB120', 'Marker', '.');
plot(snrSig, log10(varLbFreq), 'LineWidth', 2, 'Color', '#A2142F');
hold off
xlabel("SNR (dB)", "Interpreter", "latex");
ylabel("$\log_{10}(MSE_{frequency})$", "Interpreter", "latex");
ylim([-16 2]);
legend('Joint Estimator', 'Improved Joint', 'Peak Search', 'Phase Difference');
set(gca, 'Fontsize', 20);
% Plot phase MSE-SNR curve
subplot(2, 1, 2);
hold on
plot(snrSig, log10(phaMseA), 'LineWidth', 2, 'Color', '#4DBEEE', 'Marker', '*');
plot(snrSig, log10(phaMseB), 'LineWidth', 2, 'Color', '#77AC30', 'Marker', 'o');
plot(snrSig, log10(phaMseC), 'LineWidth', 2, 'Color', '#D95319', 'Marker', '+');
plot(snrSig, log10(phaMseD), 'LineWidth', 2, 'Color', '#EDB120', 'Marker', '.');
plot(snrSig, log10(varLbPha), 'LineWidth', 2, 'Color', '#A2142F');
hold off
xlabel("SNR (dB)", "Interpreter", "latex");
ylabel("$\log_{10}(MSE_{phase})$", "Interpreter", "latex");
ylim([-12 4]);
legend('Joint Estimator', 'Improved Joint', 'Peak Search', 'Phase Difference');
set(gca, 'Fontsize', 20);


%% Print Estimation Information

fprintf('\n');
fprintf('Signal frequency: %.3f Hz\n', ft);
fprintf('Signal phase: %.3f rad\n', pt);


