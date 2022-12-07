% Description:  Comparison of Estimators with Varying Sampling Time
% Projet:       Joint Estimatior of Frequency and Phase
% Date:         Nov 3, 2022
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

cycles = 0.5 : 0.2 : 2.5;         % Number of cycles
Tt = cycles / ft;                 % Total time of sampling (s)
numCycle = length(cycles);        % Iteration times
freqMseA = zeros(1, numCycle);    % MSE of frequency (Joint)
phaMseA = zeros(1, numCycle);     % MSE of phase (Joint)
timeMeanA = zeros(1, numCycle);   % Mean of time (Joint)
timeVarA = zeros(1, numCycle);    % Variance of time (Joint)
freqMseB = zeros(1, numCycle);    % MSE of frequency (Improved)
phaMseB = zeros(1, numCycle);     % MSE of phase (Improved)
timeMeanB = zeros(1, numCycle);   % Mean of time (Improved)
timeVarB = zeros(1, numCycle);    % Variance of time (Improved)
freqMseC = zeros(1, numCycle);    % MSE of frequency (Phase)
phaMseC = zeros(1, numCycle);     % MSE of phase (Phase)
timeMeanC = zeros(1, numCycle);   % Mean of time (Phase)
timeVarC = zeros(1, numCycle);    % Variance of time (Phase)
freqMseD = zeros(1, numCycle);    % MSE of frequency (Peak)
phaMseD = zeros(1, numCycle);     % MSE of phase (Peak)
timeMeanD = zeros(1, numCycle);   % Mean of time (Peak)
timeVarD = zeros(1, numCycle);    % Variance of time (Peak)

poolobj = parpool(12);
parfor i = 1 : numCycle
    
    Ns = round(Tt(i)*Fs);               % Total sampling points
    
    % Generate original signal sequence
    xt = (0 : Ns-1) / Fs;               % Time index
    at = 1;                             % Signal amplitude
    xn0 = at * sin(2*pi*ft*xt + pt);    % Test signal
    xn1 = at * cos(2*pi*ft*xt + pt);    % Test signal

    % Define estimator options
    maxIter = 5;                            % Maximum iteration time for each estimation
    numEst = 200;                           % Estimation times for each test

    % Add noise with varying SNR and estimate
    if ~noiseFlag
        xn2 = xn0;
        xn3 = xn1;
    else
        sigmaN = at / 10.^(snrSig/20);      % Standard variance of noise
        sigNoise = sigmaN * randn(1, Ns);   % Additive white Gaussian noise
        xn2 = xn0 + sigNoise;
        xn3 = xn1 + sigNoise;
    end

    % Estimate with Joint Estimator
    [freqMseA(i), phaMseA(i), timeMeanA(i), timeVarA(i)] = JointEstimatorTest(xn2, ...
        ft, pt, Fs, Tt(i), numEst, maxIter);
    % Estimate with Improved Joint Estimator
    [freqMseB(i), phaMseB(i), timeMeanB(i), timeVarB(i)] = JointEstimatorTest2(xn2, ...
        ft, pt, Fs, Tt(i), numEst, maxIter)
    % Estimate with Phase Difference Method
    [freqMseC(i), phaMseC(i), timeMeanC(i), timeVarC(i)] = PhaseDiffTest(xn3, ...
        ft, pt, Fs, Tt(i), numEst)
    % Estimate with DTFT Peak Search
    [freqMseD(i), phaMseD(i), timeMeanD(i), timeVarD(i)] = PeakSearchTest(xn3, ...
        ft, pt, Fs, Tt(i), numEst)

    fprintf('Estimation No.%d, Number of cycles = %.1f\n', i, cycles(i));

end
delete(poolobj);


%% Plot

% Plot relationship between MSE and SNR
errPlt = figure(1);
errPlt.Name = "Relationship between MSE and SNR";
errPlt.WindowState = 'maximized';
% Plot frequency MSE-SNR curve
subplot(2, 1, 1);
hold on
plot(cycles, log10(freqMseA), 'LineWidth', 2, 'Color', '#4DBEEE', 'Marker', '*');
plot(cycles, log10(freqMseB), 'LineWidth', 2, 'Color', '#77AC30', 'Marker', 'o');
plot(cycles, log10(freqMseC), 'LineWidth', 2, 'Color', '#D95319', 'Marker', '+');
plot(cycles, log10(freqMseD), 'LineWidth', 2, 'Color', '#EDB120', 'Marker', '.');
hold off
xlabel("Number of Cycles", "Interpreter", "latex");
ylabel("$\log_{10}(MSE_{frequency})$", "Interpreter", "latex");
legend('Joint Estimator', 'Improved Joint', 'Peak Search', 'Phase Difference');
set(gca, 'Fontsize', 20);
% Plot phase MSE-SNR curve
subplot(2, 1, 2);
hold on
plot(cycles, log10(phaMseA), 'LineWidth', 2, 'Color', '#4DBEEE', 'Marker', '*');
plot(cycles, log10(phaMseB), 'LineWidth', 2, 'Color', '#77AC30', 'Marker', 'o');
plot(cycles, log10(phaMseC), 'LineWidth', 2, 'Color', '#D95319', 'Marker', '+');
plot(cycles, log10(phaMseD), 'LineWidth', 2, 'Color', '#EDB120', 'Marker', '.');
hold off
xlabel("Number of Cycles", "Interpreter", "latex");
ylabel("$\log_{10}(MSE_{phase})$", "Interpreter", "latex");
legend('Joint Estimator', 'Improved Joint', 'Peak Search', 'Phase Difference');
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


