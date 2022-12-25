% Description   : Test Program for Correlation Method (AWGN, DFT-Based)
% Projet        : Joint Estimatior of Frequency and Phase
% Author        : Zhiyu Shen @Nanjing University
% Establish Date: Dec 8, 2022
% Revise Date   : Dec 8, 2022

close all
clear
clc

% Variable definition
fLb = 0;
fUb = 1;
pLb = 0;
pUb = 2*pi;
at = 1;             % Original signal's amplitude (V)
ft = fLb + 0.01*randi([0 round(100*(fUb-fLb))]);
pt = pLb + 0.01*randi([0 round(100*(pUb-pLb))]);
Fs = 5;             % Sampling rate (Hz)
Ns = 64;          % Number of sample points
Tt = 5;             % Sampling time (s)

% Generate a sine wave for test
tIdx = (0:Ns-1)/Fs;                             % Time index
sigTest = at*cos(2*pi*ft*tIdx+pt);              % Test signal

% Add gaussian white noise
SNRdB = 60;
SNRamp = 10.^(SNRdB./20);
sigmaN = at./(sqrt(2)*SNRamp);
sigNois = sigmaN*randn(1,Ns);
xn = sigTest + sigNois;

% Add window and perform DFT on the test signal
nFFT = 2^nextpow2(Ns);                          % Number of FFT points
idxWin = 0 : 1 : nFFT-1;
winSig = 0.54 - 0.46*cos(2*pi*idxWin/nFFT);     % Window signal (Hamming Window)
xnWin = [xn, zeros(1,nFFT-Ns)].*winSig;         % Zero padding and add window
xnFFT = fft(xnWin,nFFT);

% Compute the frequency spectrum of test signal
Xn1 = abs(xnFFT/nFFT);
Xn = Xn1(1:nFFT/2);
Xn(2:end-1) = 2*Xn(2:end-1);

% Sweep freqeuncy and initial phase setting
fHead = 0.01;                               % Starting frequency (Hz)
fEnd = 1;                                   % Ending frequency (Hz)
fInc = 0.0001;                              % Frequency increment (Hz)
pHead = 0;                                  % Starting phase (rad)
pEnd = 2*pi;                                % Ending phase (rad)
pInc = pi/200;                              % Phase increment (rad)

[R, iMax, jMax] = CorrSweepFreq(fHead, fEnd, fInc, pHead, pEnd, ...
    pInc, Xn, Ns, nFFT, Fs);

% Plot correlation coefficient figure
fn = 0 : iMax-1;                            % Index matrix of frequency
pn = 0 : jMax-1;                            % Index matrix of phase
idxF = fHead + fn*fInc;                     % X coordinate: frequency axis
idxP = pHead + pn*pInc;                     % Y coordinate: phase axis


%% Plot

corrPlt = figure(2);
corrPlt.Name = 'Correlation Coeffient Distribution';
corrPlt.WindowState = 'maximized';
s = surf(idxF, idxP, R.');
s.FaceAlpha = 1;
s.EdgeColor = 'flat';
s.Marker = 'none';
textParam = ['$F_s=', num2str(Fs), '\ Hz,\quad N_{sample}=', num2str(Ns), ...
    ',\quad SNR=', num2str(SNRdB), '\ dB,,\quad f_{test}=', num2str(ft), ...
    '\ Hz,\quad \phi_{test}=', num2str(pt), '\ rad$'];
title('\bf Correlation Coefficient between Measured Sequence and Constructed Sequence', ...
    textParam, 'Interpreter', 'latex', 'FontSize', 24);
xlabel('$f\ (Hz)$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\phi\ (rad)$', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('$R$', 'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'Fontsize', 20);
grid on;

% Find maximal value of correlation coefficient
[m1, m1Idx] = min(R);
[m2, m2Idx] = min(m1);
iHat = m1Idx(m2Idx)-1;
jHat = m2Idx-1;
maxR = m2;
fHat = fHead + iHat*fInc;
pHat = pHead + jHat*pInc;
fprintf('True frequency: %.4f Hz\n', ft);
fprintf('True phase: %.4f (%.4fpi) rad\n', pt, pt/pi);
fprintf('Maximal correlation coefficient: %.4f\n', maxR);
fprintf('Estimated frequency: %.4f Hz\n', fHat);
fprintf('Estimated phase: %.4f (%.4fpi) rad\n', pHat, pHat/pi);
