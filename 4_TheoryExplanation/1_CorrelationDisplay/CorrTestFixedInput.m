% Description:  Test Program for Correlation Method
%               Parameters preset
% Projet:       Joint Estimatior of Frequency and Phase
% Author:       Zhiyu Shen @Nanjing University
% Date  :       May 3, 2022

clear
clc
close all


%%% Define test signal parameters

% Variable definition
at = 1;                     % Original signal's amplitude (V)
ft = 0.02;                  % Original signal's frequency (Hz)
pt = 1/7*pi;                % Original signal's initial phase (rad)
Fs = 50;                    % Sampling rate (Hz)
Tt = 0.2 / ft;              % Sampling time (s)
addNoise = 0;

tic


%%% Generate test signal and preparation

% Generate a sine wave for test
[xnTemp, L, idx] = WaveGen(at, ft, pt, Fs, Tt);

% Add gaussian white noise
if addNoise == 0
    xn = xnTemp;
else
    SNR = 15;
    xn = awgn(xnTemp, SNR, 'measured');
end

%%% Plot

tn = idx / Fs;

% Plot time domain signal
wavePlt = figure(1);
wavePlt.Name = 'Time Domain Wave of Input Sinusoid';
wavePlt.WindowState = 'maximized';
plot(tn, xn, 'LineWidth', 2, 'Color', '#0072BD');
ylim([-1.2 * at 1.2 * at]);
textParam = ['$F_s=', num2str(Fs), '\ Hz,\quad T_{sample}=', num2str(Tt), ...
    '\ s,\quad f_{test}=', num2str(ft), '\ Hz,\quad \phi_{test}=', ...
    num2str(pt), '\ rad$'];
title('\bf Sample of Signal to Be Estimated', textParam, 'Interpreter', ...
    'latex', 'FontSize', 24);
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('Amplitude (V)', 'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'Fontsize', 20);
grid on;


%% Estimation process
% Frequency scale: 0.01 ~ 1 Hz
% Frequency precision: 0.01 Hz
% Phase scale: -pi ~ pi rad
% Phase precision: pi/100 rad
estTime = 1;
% Sweep freqeuncy and initial phase setting
a0 = 1;                                     % Amplitude (V)
fHead = 0.01;                               % Starting frequency (Hz)
fEnd = 1;                                   % Ending frequency (Hz)
fInc = 0.01;                                % Frequency increment (Hz)
pHead = -pi;                                % Starting phase (rad)
pEnd = pi;                                  % Ending phase (rad)
pInc = pi / 99;                             % Phase increment (rad)

[R, iMax, jMax] = CorrSweepTime(a0, fHead, fEnd, fInc, pHead, pEnd, pInc, ...
    idx, xn, Fs);

% Plot correlation coefficient figure
fn = 0 : iMax;                              % Index matrix of frequency
pn = 0 : jMax;                              % Index matrix of phase
fIdx = fHead + fn * fInc;                   % X coordinate: frequency axis
pIdx = pHead + pn * pInc;                   % Y coordinate: phase axis
[X, Y] = meshgrid(fIdx, pIdx);              % Generate coordinate grid
Z = (R(fn + 1, pn + 1)).';                  % Z coordinate: cross correlation coefficient axis

corrPlt = figure(2);
corrPlt.Name = 'Correlation Coeffient Distribution';
corrPlt.WindowState = 'maximized';
s = surf(X, Y, Z);
s.FaceAlpha = 1;
s.EdgeColor = 'flat';
s.Marker = 'none';
textParam = ['$F_s=', num2str(Fs), '\ Hz,\quad T_{sample}=', num2str(Tt), ...
    '\ s,\quad f_{test}=', num2str(ft), '\ Hz,\quad \phi_{test}=', ...
    num2str(pt), '\ rad$'];
title('\bf Correlation Coefficient between Measured Sequence and Constructed Sequence', ...
    textParam, 'Interpreter', 'latex', 'FontSize', 24);
xlabel('$f\ (Hz)$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\phi\ (rad)$', 'Interpreter', 'latex', 'FontSize', 20);
zlabel('$R$', 'Interpreter', 'latex', 'FontSize', 20);
set(gca, 'Fontsize', 20);
grid on;

