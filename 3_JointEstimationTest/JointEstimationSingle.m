% Description:  Test Program for Joint Estimator for Single Run
% Projet:       Short Sequence Parameter Estimation
% Author:       Zhiyu Shen @Nanjing University
% Date  :       July 28, 2022

clear
close all
clc

Fs = 100;                           % Sampling frequency (Hz)
Tt = 20;                            % Total time of sampling (s)
Ns = Tt * Fs;                       % Total sampling points

% ft = randi([8 100]) / 100;              % Frequency of test signal (Hz)
% pt = (randi([0 200]) - 100) * pi / 100; % Phase of test signal (rad)

ft = 0.2;                               % Frequency of test signal (Hz)
pt = pi / 3;                            % Phase of test signal (rad)
at = 1;                                 % Amplitude of test signal (V)
xSig = [ft; pt; at];

xt = (0 : Ns - 1) / Fs;                 % Time index
xn = at * sin(2 * pi * ft * xt + pt);   % Test signal

M = 50;                             % Search times

options.maxIter = M;
tic
[xBest, yBest, info] = JointEstimator(xn, Fs, options);
toc

estErr = abs(xBest - xSig) ./ xSig;
fErr = estErr(1);
pErr = estErr(2);
aErr = estErr(3);

fprintf('\n-------- Input Signal --------\n');
fprintf('Frequency: %.3d Hz\n', ft);
fprintf('Phase: %.3d rad\n', pt);

fprintf('\n-------- Estimation Result --------\n');
fprintf('Frequency: %.3d Hz\n', fe);
fprintf('Phase: %.3d rad\n', pe);

fprintf('\n-------- Error Analysis --------\n');
fprintf('Frequency Error: %.3d\n', fErr);
fprintf('Phase Error: %.3d\n', pErr);


