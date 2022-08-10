% Description:  Test Program for Joint Estimator for Changing Frequence
% Projet:       Short Sequence Parameter Estimation
% Date:         Aug 6, 2022
% Author:       Zhiyu Shen

clear
close all
clc

Fs = 100;                           % Sampling frequency (Hz)
Tt = 2;                             % Total time of sampling (s)
Ns = Tt * Fs;                       % Total sampling points

df = 0.01;                          % Frequency increment
fmin = 0.05;                        % Minimum frequency
fmax = 1;                           % Maximum frequency
imax = (fmax - fmin) / df + 1;      % Maximum loop index

fIdx = 1 : imax;
ft = fmin + (fIdx - 1) * df;
pt = (randi([1 19]) - 10) * pi / 10;    % Phase of test signal (rad)

M = 50;                             % Search times
options.maxIter = M;

fErr = zeros(1, imax);
pErr = zeros(1, imax);
iterTime = zeros(1, imax);

p = parpool(8);
parfor i = 1 : imax

    xt = (0 : Ns - 1) / Fs;                 % Time index
    xn = sin(2 * pi * ft(i) * xt + pt);     % Test signal

    tic
    [xBest, ~, ~] = JointEstimator(xn, Fs, options);
    iterTime(i) = toc;

    fe = xBest(1);
    pe = xBest(2);
    fErr(i) = abs(fe - ft(i)) / ft(i);
    pErr(i) = abs(pe - pt) / pt;

end

delete(p);

totTime = sum(iterTime);
abgTime = totTime / imax;

fprintf('\nCompleted!\n');
fprintf('Total Run Time: %.2d\n', totTime);
fprintf('\n-------- Input Signal --------\n');
fprintf('Frequency range: %.3d ~ %d Hz\n', fmin, fmax);
fprintf('Phase: %.3d rad\n', pt);

figure(1)
semilogy(ft, fErr, "LineWidth", 2, "Color", "#0072BD", "Marker", "square");



