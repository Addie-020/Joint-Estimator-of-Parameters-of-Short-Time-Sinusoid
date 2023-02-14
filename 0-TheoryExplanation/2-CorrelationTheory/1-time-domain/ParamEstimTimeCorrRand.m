% Description   :  Test Program for Correlation Method (Time Domain, AWGN)
% Projet        :  Joint Estimatior of Frequency and Phase
% Author        :  Zhiyu Shen @Nanjing University
% Establish Date:  Dec 4, 2022
% Revise Date   :  Dec 12, 2022

close all
clear
clc

% Variable definition
fLb = 0.5;
fUb = 1.5;
pLb = 0;
pUb = 2*pi;
at = 1;                 % Original signal's amplitude (V)
ft = fLb + 0.01*randi([0 round(100*(fUb-fLb))]);
pt = pLb + 0.01*randi([0 round(100*(pUb-pLb))]);
Fs = 50;                % Sampling rate (Hz)
Ns = 64;                % Number of samples

% Generate a sine wave for test
tIdx = (0:Ns-1)/Fs;                             % Time index
sigTest = at*cos(2*pi*ft*tIdx+pt);              % Test signal

% Add gaussian white noise
SNRdB = 8;
snr = 10.^(SNRdB./20);
sigmaN = at./(sqrt(2)*snr);
sigNois = sigmaN*randn(1,Ns);
sigMeas = sigTest + sigNois;

% Sweep freqeuncy and initial phase setting
fHead = fLb;                               % Starting frequency (Hz)
fEnd = fUb;                                   % Ending frequency (Hz)
fInc = 0.001;                              % Frequency increment (Hz)
pHead = pLb;                                  % Starting phase (rad)
pEnd = pUb;                                % Ending phase (rad)
pInc = pi/100;                              % Phase increment (rad)

[R, iMax, jMax] = CorrSweepTime(fHead, fEnd, fInc, pHead, pEnd, ...
    pInc, Ns, sigMeas, Fs);

% Plot correlation coefficient figure
fn = 0 : iMax-1;                            % Index matrix of frequency
pn = 0 : jMax-1;                            % Index matrix of phase
idxF = fHead + fn*fInc;                     % X coordinate: frequency axis
idxP = pHead + pn*pInc;                     % Y coordinate: phase axis

corrPlt = figure(2);
corrPlt.Name = 'Correlation Coeffient Distribution';
corrPlt.WindowState = 'maximized';
s = surf(idxF, idxP/pi, R.');
s.FaceAlpha = 1;
s.EdgeColor = 'flat';
s.Marker = 'none';
textParam = ['$F_s=', num2str(Fs), '\ Hz,\quad N_{sample}=', num2str(Ns), ...
    ',\quad SNR=', num2str(SNRdB), '\ dB,,\quad f_{test}=', num2str(ft), ...
    '\ Hz,\quad \phi_{test}=', num2str(pt), '\ rad$'];
title('\bf Correlation Coefficient between Measured Sequence and Constructed Sequence', ...
    textParam, 'Interpreter', 'latex', 'FontSize', 24);
xlabel('$f\ (Hz)$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\phi\ (\pi rad)$', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('$R$', 'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'Fontsize', 20);
grid on;
colorbar

% Find maximal value of correlation coefficient
[m1, m1Idx] = max(R);
[m2, m2Idx] = max(m1);
iHat = m1Idx(m2Idx)-1;
jHat = m2Idx-1;
maxR = m2;
fHat = fHead + iHat*fInc;
pHat = pHead + jHat*pInc;
fprintf('SNR = %.2f dB = %.2f\n', SNRdB, snr);
fprintf('True frequency: %.4f Hz\n', ft);
fprintf('True phase: %.4f (%.4fpi) rad\n', pt, pt/pi);
fprintf('Maximal correlation coefficient: %.4f\n', maxR);
fprintf('Estimated frequency: %.4f Hz\n', fHat);
fprintf('Estimated phase: %.4f (%.4fpi) rad\n', pHat, pHat/pi);
