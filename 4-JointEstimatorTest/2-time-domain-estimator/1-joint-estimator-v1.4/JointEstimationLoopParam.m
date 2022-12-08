% Description:  Test Program for Joint Estimator for Varying Signal Parameters
% Projet:       Joint Estimatior of Frequency and Phase
% Date:         Nov 3, 2022
% Author:       Zhiyu Shen

clear
close all
clc

%% Set Noise Figure

addNoise = input('Add noise to signal? Y/N [N]: ', 's');
if isempty(addNoise) || (addNoise == 'N')
    noiseFlag = 0;
    snrSig = 1000;
elseif addNoise == 'Y'
    % Define SNR
    snrSig = input('SNR(dB) [40]: ');
    if isempty(snrSig)
        snrSig = 40;
    end
    noiseFlag = 1;
end

cycles = input('Number of cycles sampled: [0.5]: ');
if isempty(cycles)
    cycles = 0.5;
end


%% Iteration

% Define parameters
ft = 0.01 : 0.01 : 0.5;                 % Frequency range
pt = 0 : pi/50 : pi*99/50;              % Phase range
at = 1;                                 % Signal amplitude
numFreq = length(ft);
numPha = length(pt);

% Generate test signals
Tt = cycles ./ ft;                      % Total time of sampling (s)
Fs = 10;                                % Sampling frequency (Hz)

% Allocate memory for iteration parameters
freqMse = zeros(numFreq, numPha);       % MSE of frequency
phaMse = zeros(numFreq, numPha);        % MSE of phase
timeMean = zeros(numFreq, numPha);      % Mean of time
timeVar = zeros(numFreq, numPha);       % Variance of time

% Define estimator options
options.maxIter = 5;                    % Maximum iteration time for each estimation
numEst = 1;

% Outer loop: frequency
for i = 1 : numFreq
    
    % Generate original signal sequence parameters
    Ns = round(Tt(i)*Fs);               % Total sampling points
    xt = (0 : Ns-1) / Fs;               % Time index
    
    % Add noise with varying SNR and estimate
    if ~noiseFlag
        sigNoise = zeros(1, Ns);
    else
        sigmaN = at / 10.^(snrSig/20);      % Standard variance of noise
        sigNoise = sigmaN * randn(1, Ns);   % Additive white Gaussian noise
    end
    
    % Inner loop: phase
    for j = 1: numPha
        xn0 = at*sin(2*pi*ft(i)*xt+pt(j));  % Original signal
        xn = xn0 + sigNoise;                % Add noise    

        [freqMse(i,j), phaMse(i,j), timeMean(i,j), timeVar(i,j)] = ...
            JointEstimatorTest(xn, ft(i), pt(j), Fs, Tt(i), numEst, options, [], []);
    end

    fprintf('Estimation No.%d, Frequency = %.2f Hz\n', i, ft(i));

end


%% Plot

% Text displayed
textParam = ['$F_s=', num2str(Fs), '$ Hz, $t_s=', num2str(cycles), ...
    'T$, $SNR=', num2str(snrSig), '$ dB'];

% Plot variance of mean estimation time
timePlt = figure(1);
timePlt.Name = "Variance of Mean Estimation Time";
timePlt.WindowState = 'maximized';
% Plot curve
surfTime = surf(pt, ft, timeMean);
surfTime.FaceAlpha = 0.8;
surfTime.EdgeColor = 'interp';
surfTime.Marker = 'none';
colorbar
% Set the plotting properties
xlabel("Phase $\phi\ (rad)$", "Interpreter", "latex");
ylabel("Frequency $f\ (Hz)$", "Interpreter", "latex");
zlabel("Mean Estimation Time (s)", "Interpreter", "latex");
set(gca, 'Fontsize', 20);

% Plot variance of frequency MSE
fmsePlt = figure(2);
fmsePlt.Name = "Variance of Frequency MSE";
fmsePlt.WindowState = 'maximized';
% Plot curve
surfFreq = surf(pt, ft, log10(freqMse));
surfFreq.FaceAlpha = 0.8;
surfFreq.EdgeColor = 'interp';
surfFreq.Marker = 'none';
colorbar
% Set the plotting properties
xlabel("Phase $\phi\ (rad)$", "Interpreter", "latex")
ylabel("Frequency $f\ (Hz)$", "Interpreter", "latex");
zlabel("$\log_{10}(MSE_{frequency})$", "Interpreter", "latex");
xlim([0 max(pt)]);
ylim([0 max(ft)]);
set(gca, 'Fontsize', 20);

% Plot variance of phase MSE
fmsePlt = figure(3);
fmsePlt.Name = "Variance of Phase MSE";
fmsePlt.WindowState = 'maximized';
% Plot curve
surfPha = surf(pt, ft, log10(phaMse));
surfPha.FaceAlpha = 0.8;
surfPha.EdgeColor = 'interp';
surfPha.Marker = 'none';
colorbar
% Set the plotting properties
xlabel("Phase $\phi\ (rad)$", "Interpreter", "latex")
ylabel("Frequency $f\ (Hz)$", "Interpreter", "latex");
zlabel("$\log_{10}(MSE_{phase})$", "Interpreter", "latex");
xlim([0 max(pt)]);
ylim([0 max(ft)]);
set(gca, 'Fontsize', 20);


%% Print Estimation Information

if ~noiseFlag
    fprintf('Noise not added.\n');
else
    fprintf('SNR: %.3f dB\n', snrSig);
end



