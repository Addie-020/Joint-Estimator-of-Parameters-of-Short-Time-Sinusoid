% Description:  Test program for joint estimator for changing arameters
%               Both frequency and initial phase are changing
% Projet:       Short Sequence Parameter Estimation
% Date:         Aug 22, 2022
% Author:       Zhiyu Shen

clear
close all
clc

%%% Parameters Definition

% Basic parameters
Fs = 100;                           % Sampling frequency (Hz)

% Frequency changing parameters
df = 0.01;                          % Frequency increment
fmin = 0.01;                        % Minimum frequency
fmax = 0.1;                         % Maximum frequency
imax = fix((fmax - fmin) / df + 1); % Maximum frequency loop index

% Phase changing parameters
dp = pi/19;                         % Phase increment
pmin = pi/10;                       % Minimum phase
pmax = 2*pi;                        % Maximum phase
jmax = fix((pmax - pmin) / dp + 1); % Maximum phase loop index


%%% Estimation Process

% Create parameter vectors
fIdx = 1 : imax;
pIdx = 1 : jmax;
ft = fmin + (fIdx - 1) * df;        % Frequency of test signal (Hz)
pt = pmin + (pIdx - 1) * dp;        % Phase of test signal (rad)

% Set search times
M = 50;                             % Search times
options.maxIter = M;

% Allocate memory for vectors storing estimation results
fErr = zeros(imax, jmax);
pErr = zeros(imax, jmax);
iterTime = zeros(imax, jmax);

% Main loop
for i = 1 : imax
    for j = 1 : jmax

        Tt = 0.1 / ft(i);                       % Total time of sampling (s)
        Ns = Tt * Fs;                           % Total sampling points

        xt = (0 : Ns - 1) / Fs;                 % Time index
        xn = sin(2 * pi * ft(i) * xt + pt(j));  % Test signal

        tic
        [xBest, ~, ~] = JointEstimator(xn, Fs, options);
        iterTime(i, j) = toc;

        fe = xBest(1);
        pe = xBest(2);
        fErr(i, j) = abs((fe - ft(i)) / ft(i));
        pErr(i, j) = abs((pe - pt(j)) / pt(j));

    end
end


%%% Display the Estimation Result

% Calculate computation time
totTime = sum(iterTime, "all");
avgTime = totTime / (imax * jmax);

% Display estimation infomation
fprintf('\nEstimation Completed!\n');
fprintf('Total Run Time: %.2d\n', totTime);
fprintf('Average Run Time: %.2d\n', avgTime);
fprintf('\n-------- Input Signal --------\n');
fprintf('Frequency range: %.3d ~ %.3d Hz\n', fmin, fmax);
fprintf('Phase range: %.3d ~ %.3d rad\n', pmin, pmax);


%%% Plot Estimation Error

% Figure settings
errPlt = figure(1);
errPlt.WindowState = 'maximized';

% Calculate data index
[X, Y] = meshgrid(ft, pt);              % Generate coordinate grid
Z1 = fErr(fIdx, pIdx).';                % Z coordinate: frequency estimation error axis
Z2 = pErr(fIdx, pIdx).';                % Z coordinate: phase estimation error axis

% Plot frequency error
subplot(2, 1, 1);
s1 = surf(X, Y, Z1);
s1.FaceAlpha = 1;
s1.EdgeColor = 'flat';
s1.Marker = 'none';
title('Frequency Estimation Error');
xlabel('\it f /\rm Hz');
ylabel('\it \phi /\rm rad');
zlabel('frequency error');
zlim([0 1]);
set(gca, 'FontSize', 10 , 'FontName', 'Times New Roman');
grid on;

% Plot phase error
subplot(2, 1, 2);
s2 = surf(X, Y, Z2);
s2.FaceAlpha = 1;
s2.EdgeColor = 'flat';
s2.Marker = 'none';
title('Phase Estimation Error');
xlabel('\it f /\rm Hz');
ylabel('\it \phi /\rm rad');
zlabel('phase error');
zlim([0 1]);
set(gca, 'FontSize', 10 , 'FontName', 'Times New Roman');
grid on;



