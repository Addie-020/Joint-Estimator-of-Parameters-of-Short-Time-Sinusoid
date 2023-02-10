% Description:  Test Program for Verifying Theory of Power Spectrum
% Projet:       Joint Estimatior of Frequency and Phase
% Date:         Feb 1, 2023
% Author:       Zhiyu Shen

clear
close all
clc


%% Test Signal Construction

% Define signal parameters
numSig = 3;                 % Number of test signals
a0 = 1;                     % Amplitude of test signal (V)
f0 = [0.5 0.5 0.5];                 % Frequency of test signal (Hz)
p0 = [pi/8 0 pi*1/2];              % Phase of test signal (rad)

% Define system parameters
Fs = 5;                     % Sampling rate (Hz)
Ts = 1/Fs;                  % Sampling interval (s)
numSample = 32;             % Number of sampling points

% Generate discrete-time signal
idxT = (0:numSample-1).'*Ts;                        % Time index of sampling points
xn = a0*cos(2*pi*idxT*f0+repmat(p0,numSample,1));   % Discrete-time test signal


%% Original Signal PSD Computation

% Add window to signal
win = boxcar(numSample);                % Generate a window signal
xnWin = xn.*repmat(win,1,numSig);       % Add window to test signal

% Compute FFT sequence of test signal
numFFT = 2^nextpow2(numSample);         % Calculate number of FFT ponits
xdft = fft(xnWin,numFFT,1);             % Compute FFT of windowed test signal
xdft = xdft(1:numFFT/2+1,:);            % Retain the values of positive frequency

% Compute PSD based on FFT values
psdx = 1/(Fs*numFFT) * abs(xdft).^2;    % Compute PSD of test signal
psdx(2:end-1,:) = psdx(2:end-1,:)*2;    % Normalize the PSD
idxF = 0:Fs/numFFT:Fs/2;                % Compute frequency index of PSD


%% Add noise to signal

% Define signal parameters (White Gaussian noise)
SNRdB = 10;                             % SNR of power in dB
sigma = sqrt(a0.^2/10.^(SNRdB/10));     % Standard variance of noise
wn = sigma*randn(numSample,1);          % Gaussian noise sequence

% Add noise to signal
yn = xn+repmat(wn,1,numSig);

% Add window
ynWin = yn.*repmat(win,1,numSig);

% Compute PSD
ydft = fft(ynWin,numFFT,1);
ydft = ydft(1:numFFT/2+1,:);
psdy = 1/(Fs*numFFT) * abs(ydft).^2;
psdy(2:end-1,:) = psdy(2:end-1,:)*2;


%% Compute PSD of noise signal

% Add window
wnWin = wn.*repmat(win,1,numSig);

% Compute PSD
wdft = fft(wnWin,numFFT,1);
wdft = wdft(1:numFFT/2+1,:);
psdw = 1/(Fs*numFFT) * abs(wdft).^2;
psdw(2:end-1,:) = psdw(2:end-1,:)*2;



%% Plot

figure

% Plot PSD of original signals
subplot(3,1,1);
hold on
plot(idxF,psdx(:,1), 'LineWidth', 1, 'Color', '#0072BD');
plot(idxF,psdx(:,2), 'LineWidth', 1, 'Color', '#D95319');
plot(idxF,psdx(:,3), 'LineWidth', 1, 'Color', '#EDB120');
legend('$0.1$ Hz, $\pi/3$ rad','$0.1$ Hz, $\pi/2$ rad', ...
    '$0.1$ Hz, $2\pi/5$ rad','Interpreter','latex');
hold off

% Plot PSD of noise signal
subplot(3,1,2);
hold on
plot(idxF,psdw,'LineWidth',1,'Color','#0072BD');
hold off

% Plot PSD of signal with noise
subplot(3,1,3);
hold on
plot(idxF,psdy(:,1),'LineWidth',1,'Color','#0072BD');
plot(idxF,psdy(:,2),'LineWidth',1,'Color','#D95319');
plot(idxF,psdy(:,3),'LineWidth',1,'Color','#EDB120');
legend('$0.1$ Hz, $\pi/3$ rad','$0.1$ Hz, $\pi/2$ rad', ...
    '$0.1$ Hz, $2\pi/5$ rad','Interpreter','latex');
hold off


%% Compute Correlation Coefficients

objFunMat = zeros(numSig,numSig);
for i = 1:numSig
    for j = 1:numSig
        corrTemp = corrcoef(psdx(:,i),psdy(:,j));
        corrVal = corrTemp(1,2);
        objFunValTemp = exp((1-corrVal)*1e3);
        objFunMat(i,j) = min(objFunValTemp,1e5);
    end
end