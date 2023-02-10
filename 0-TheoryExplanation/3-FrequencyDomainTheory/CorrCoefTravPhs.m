% Description   : Visualization of Correlation Coefficient with Varying Phase
% Projet        : Joint Estimatior of Frequency and Phase
% Author        : Zhiyu Shen @Nanjing University
% Establish Date: Dec 8, 2022
% Revise Date   : Feb 3, 2022

close all
clear
clc


%% Test Signal Construction

% Define signal parameters
a0 = 1;                     % Amplitude of test signal (V)
f0 = 0.1;                   % Frequency of test signal (Hz)
p0 = pi/3;                    % Phase of test signal (rad)

% Define system parameters
Fs = 5;                     % Sampling rate (Hz)
Ts = 1/Fs;                  % Sampling interval (s)
numSample = 128;            % Number of sampling points

% Generate discrete-time signal
idxT = (0:numSample-1).'*Ts;                        % Time index of sampling points
xn = a0*cos(2*pi*idxT*f0+p0);   % Discrete-time test signal


%% Original Signal PSD Computation

% Add window to signal
% win = boxcar(numSample);                % Generate a window signal
win = hamming(numSample);
xnWin = xn.*win;                        % Add window to test signal

% Compute FFT sequence of test signal
numFFT = 2^nextpow2(numSample);         % Calculate number of FFT ponits
xdft = fft(xnWin,numFFT,1);             % Compute FFT of windowed test signal
xdft = xdft(1:numFFT/2+1,:);            % Retain the values of positive frequency

% Compute PSD based on FFT values
psdx = 1/(Fs*numFFT) * abs(xdft).^2;    % Compute PSD of test signal
psdx(2:end-1,:) = psdx(2:end-1,:)*2;    % Normalize the PSD


%% Construct Sweep Signals

% Sweep phase range
pHead = 0;                                  % Starting phase (rad)
pEnd = 2*pi;                                % Ending phase (rad)
pInc = pi/100;                              % Phase increment (rad)

% Sweep signal parameters
a1 = a0;
f1 = f0;
p1 = pHead:pInc:pEnd;
numPha = length(p1);

% Construct sweep signals
yn = a1*cos(2*pi*f1*idxT+p1);


%% Sweep Signals PSD Computation

% Add window to signal
ynWin = yn.*win;                        % Add window to test signal

% Compute FFT sequence of sweep signals
ydft = fft(ynWin,numFFT,1);             % Compute FFT of windowed test signal
ydft = ydft(1:numFFT/2+1,:);            % Retain the values of positive frequency

% Compute PSD based on FFT values
psdy = 1/(Fs*numFFT) * abs(ydft).^2;    % Compute PSD of test signal
psdy(2:end-1,:) = psdy(2:end-1,:)*2;    % Normalize the PSD


%% Compute Correlation Coefficients

% Allocate memory for vectors storing results
corrVal = zeros(1,numPha);
corrValComp = zeros(1,numPha);

% Sweep phase and compute coefficients based on power spectrum
for i = 1:numPha
    corrVal(i) = CorrCal(psdx, psdy(:,i));
end

% Sweep phase and compute coefficients based on DFT value
for i = 1:numPha
    corrValComp(i) = ComplexCorrCal(xdft, ydft(:,i));
end


%% Plot

figure

subplot(2,1,1)
hold on
plot(p1,corrVal, 'LineWidth', 1, 'Color', '#0072BD');
hold off

subplot(2,1,2)
hold on
plot(p1,corrValComp, 'LineWidth', 1, 'Color', '#0072BD');
hold off