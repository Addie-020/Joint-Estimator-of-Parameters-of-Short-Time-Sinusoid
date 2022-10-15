% Description:  Test Program for Joint Estimator for Single Run
% Projet:       Joint Estimatior of Frequency and Phase
% Date:         Sept 28, 2022
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
Tt = 0.3 / ft;                      % Total time of sampling (s)
Ns = round(Tt * Fs);                % Total sampling points

% Generate original signal sequence
xt = (0 : Ns - 1) / Fs;             % Time index
at = 1;                             % Signal amplitude
xn0 = at * sin(2*pi*ft*xt + pt);    % Test signal

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
numEst = 10;                        % Estimation times
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
timeMax = max(timeEst);
timeMin = min(timeEst);
timeMean = sum(timeEst) ./ numEst;
timeVar = sum((timeEst-timeMean).^2) / numEst;

% Calculate error
freqErr = sum((fe-ft).^2) ./ numEst;
phaErr = sum((pe-pt).^2) ./ numEst;


%%% Display and plot

% Display information
fprintf('\n-------- Input Signal --------\n');
fprintf('Frequency: %.3d Hz\n', ft);
fprintf('Phase: %.3d rad\n', pt);

fprintf('\n-------- Time Used --------\n');
fprintf('Sampling time: %.3f s\n', Tt);
fprintf('Minimum time: %.3f s\n', timeMin);
fprintf('Maximum time: %.3f s\n', timeMax);
fprintf('Mean of time: %.3f s\n', timeMean);
fprintf('Variance of time: %.5f s\n', timeVar);

fprintf('\n-------- Estimation Result --------\n');
fprintf('Frequency MSE: %.3d\n', freqErr);
fprintf('Phase MSE: %.3d\n', phaErr);

% Plot

% Plot relationship between time and iteration
timePlt = figure(1);
timePlt.Name = "Relationship between Time and Iteration";
timePlt.WindowState = 'maximized';
% Plot curve
hold on
plot(1:numEst, timeEst, 'LineWidth', 2, 'Color', '#D95319', 'Marker', '*', 'MarkerSize', 8);
hold off
% Set the plotting properties
xlabel("Iteration Index", "Interpreter", "latex");
ylabel("Estimation Time (s)", "Interpreter", "latex");
ylim([Tt-2, Tt+2]);
set(gca, 'Fontsize', 20);



