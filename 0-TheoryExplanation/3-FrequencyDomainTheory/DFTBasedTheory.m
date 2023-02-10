% Description:  Test Program for Verifying DFT Based Theory
% Projet:       Joint Estimatior of Frequency and Phase
% Date:         Feb 2, 2023
% Author:       Zhiyu Shen

clear
close all
clc


%% Test Signal Construction

% Define signal parameters
numSig = 3;                 % Number of test signals
a0 = 1;                     % Amplitude of test signal (V)
f0 = [1 1 1];                 % Frequency of test signal (Hz)
p0 = [pi/9 pi*7/4 pi*1/2-1];              % Phase of test signal (rad)

% Define system parameters
Fs = 5;                     % Sampling rate (Hz)
Ts = 1/Fs;                  % Sampling interval (s)
numSample = 64;            % Number of sampling points

% Generate discrete-time signal
idxT = (0:numSample-1).'*Ts;                        % Time index of sampling points
xn = a0*cos(2*pi*idxT*f0+repmat(p0,numSample,1));   % Discrete-time test signal


%% Original Signal PSD Computation

% Add window to signal
win = boxcar(numSample);                % Generate a window signal
xnWin = xn.*repmat(win,1,numSig);       % Add window to test signal

% Compute DFT sequence of test signal
numFFT = 2^nextpow2(numSample);         % Calculate number of FFT ponits
xdft = fft(xnWin,numFFT,1);             % Compute FFT of windowed test signal
xdft = xdft(1:numFFT/2+1,:);            % Retain the values of positive frequency


%% Add noise to signal

% Define signal parameters (White Gaussian noise)
SNRdB = 100;                             % SNR of power in dB
sigma = sqrt(a0.^2/10.^(SNRdB/10));     % Standard variance of noise
wn = sigma*randn(numSample,1);          % Gaussian noise sequence

% Add noise to signal
yn = xn+repmat(wn,1,numSig);

% Add window
ynWin = yn.*repmat(win,1,numSig);

% Compute DFT
ydft = fft(ynWin,numFFT,1);
ydft = ydft(1:numFFT/2+1,:);


%% Compute Correlation Coefficients

corrVal = zeros(numSig,numSig);
objFunMat = zeros(numSig,numSig);
for i = 1:numSig
    for j = 1:numSig
        corrVal(i,j) = CorrCal(xdft(:,i),ydft(:,j));
        objFunValTemp = exp((1-corrVal(i,j))*1e3);
        objFunMat(i,j) = min(objFunValTemp,1e5);
    end
end