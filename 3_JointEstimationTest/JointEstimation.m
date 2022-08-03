clear
close all
clc

Fs = 100;                           % Sampling frequency (Hz)
Tt = 2;                             % Total time of sampling (s)
Ns = Tt * Fs;                       % Total sampling points

ft = 0.2;                           % Frequency of test signal (Hz)
pt = pi / 3;                        % Phase of test signal (rad)

xt = (0 : Ns - 1) / Fs;             % Time index
xn = sin(2 * pi * ft * xt + pt);    % Test signal

M = 50;                             % Search times

options.maxIter = M;
tic
[xBest, yBest, info] = JointEstimator(xn, Fs, options);
toc


