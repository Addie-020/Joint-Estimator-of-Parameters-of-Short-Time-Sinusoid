% Description:  Test Program for Joint Estimator for Single Run
% Projet:       Joint Estimatior of Frequency and Phase
% Date:         July 28, 2022
% Author:       Zhiyu Shen

clear
close all
clc

% ft = randi([8 100]) / 100;              % Frequency of test signal (Hz)
% pt = (randi([0 200]) - 100) * pi / 100; % Phase of test signal (rad)

ft = 0.01;                          % Frequency of test signal (Hz)
pt = pi/3;                          % Phase of test signal (rad)

Fs = 10;                            % Sampling frequency (Hz)
Tt = 0.3 / ft;                      % Total time of sampling (s)
Ns = round(Tt * Fs);                % Total sampling points

xt = (0 : Ns - 1) / Fs;             % Time index
at = 1;                             % Signal amplitude
xn0 = at * sin(2*pi*ft*xt + pt);    % Test signal

snrSig = 40;                        % SNR (dB)
sigmaN = at / 10.^(snrSig/20);      % Standard variance of noise
sigNoise = sigmaN * randn(1, Ns);   % Additive white Gaussian noise

xn = xn0 + sigNoise;

M = 10;                             % Search times

options.maxIter = M;
tic
[xBest, yBest, info] = JointEstimator(xn, Fs, options);
totTime = toc;

fe = xBest(1);
pe = xBest(2);
fErr = abs((fe-ft) / ft);
pErr = abs((pe-pt) / pt);

fprintf('\n-------- Input Signal --------\n');
fprintf('Frequency: %.3d Hz\n', ft);
fprintf('Phase: %.3d rad\n', pt);

fprintf('\n-------- Estimation Result --------\n');
fprintf('Frequency: %.3d Hz\n', fe);
fprintf('Phase: %.3d rad\n', pe);

fprintf('\n-------- Error Analysis --------\n');
fprintf('Frequency Error: %.3d\n', fErr);
fprintf('Phase Error: %.3d\n', pErr);

fprintf('\n-------- Time Used --------\n');
fprintf('Sampling time: %.3f s\n', Tt);
fprintf('Total time: %.3f s\n', totTime);
fprintf('Mean time: %.3f s\n', info.meanTime);
fprintf('Test time: %.3f s\n', totTime+Tt);


