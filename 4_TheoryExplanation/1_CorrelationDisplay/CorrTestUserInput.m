% Description:  Test Program for Correlation Method
%               Parameters from keyboard input
% Projet:       Short Sequence Parameter Estimation
% Author:       Zhiyu Shen @Nanjing University
% Date  :       May 3, 2022

clear
clc

fprintf('-------- Short Sequence Signal Frequency and Initial Phase Co-Estimation --------\n\n');

row = 2;
line = 2;

%%% Fetch test signal parameters from keyboard input

fprintf('Please enter parameter values of test signal...\n');
at = input('Amplitude (1V): ');
ft = input('Frequency (0 ~ 1 Hz): ');
pt = input('Initial phase (-pi ~ pi rad): ');
fprintf('Please input sampling parameters...\n');
Fs = input('Sampling rate (Hz): ');
Tt = input('Sampling time (s): ');
addNoise = input('Adding noise? (1: Yes, 0: No): ');
fprintf('The test signal is %.0f*sin(2pi*%.2f+%.2f)\n', at, ft, pt);
fprintf('Signal sampling rate is %.2fHz and sample for %.2fs\n', Fs, Tt);

tic

%%% Generate test signal and preparation

fprintf('\n------------------------------- Estimation Start -------------------------------\n');

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
plot(tn, xn, 'b-', 'LineWidth', 1.5);
ylim([-1.2 * at 1.2 * at]);
title('Time domain signal');
xlabel('\it t /\rm s');
ylabel('\it f \rm(\itt\rm)');
set(gca, 'FontSize', 16);
grid on;


%% Estimation process 1
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

[R, iMax, jMax] = correlation_sweep_time(a0, fHead, fEnd, fInc, pHead, pEnd, pInc, ...
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
title('Correlation Coefficient');
xlabel('\it f /\rm Hz');
ylabel('\it \phi /\rm rad');
zlabel('\it R \rm(\itf,\phi\rm)');
set(gca, 'FontSize', 16);
grid on;

% Find maximal value of correlation coefficient
rTemp = R;
[m, idxM] = max(rTemp);
[m1, idxM1] = max(m);
iHat = idxM(idxM1) - 1;
jHat = idxM1 - 1;
rMax = m1;
fHat = fHead + iHat * fInc;
pHat = pHead + jHat * pInc;
fprintf('Estimation No.%d:\n', estTime);
fprintf('Precision: %.5f Hz, %.5f pi rad\n', fInc, pInc/pi);
fprintf('Maximal correlation coefficient: %.5f\n', rMax);
fprintf('Frequency: %.5f Hz\n', fHat);
fprintf('Phase: %.5f rad (%.5fpi)\n', pHat, pHat/pi);

%% Estimation process 2
% Frequency scale: 0.(x-1) ~ 0.(x+1) Hz
% Precision: 0.002 Hz
% Phase scale: 0.(x-1)*pi ~ 0.(x+1)*pi rad
% Phase precision: pi/500 rad
estTime = estTime + 1;
fe = round(fHat * 10);
pe = round(pHat / pi * 10);

% Sweep freqeuncy and initial phase setting
a0 = 1;                                     % Amplitude (V)
fHead = max(0, (fe - 1)) / 10;              % Starting frequency (Hz)
fEnd = (fe + 1) / 10;                       % Ending frequency (Hz)
fInc = 0.002;                               % Frequency increment (Hz)
pHead = (pe - 1) / 10 * pi;                 % Starting phase (rad)
pEnd = (pe + 1) / 10 * pi;                  % Ending phase (rad)
pInc = pi / 499;                            % Phase increment (rad)

[R, iMax, jMax] = correlation_sweep_time(a0, fHead, fEnd, fInc, pHead, pEnd, pInc, ...
    idx, xn, Fs);

% Plot correlation coefficient figure
iterPlt = figure(3);
iterPlt.Name = 'Correlation Coeffient Distribution of Each Iteration';
iterPlt.WindowState = 'maximized';
fn = 0 : iMax;                              % Index matrix of frequency
pn = 0 : jMax;                              % Index matrix of phase
fIdx = fHead + fn * fInc;                   % X coordinate: frequency axis
pIdx = pHead + pn * pInc;                   % Y coordinate: phase axis
[X, Y] = meshgrid(fIdx, pIdx);              % Generate coordinate grid
Z = (R(fn + 1, pn + 1)).';                  % Z coordinate: cross correlation coefficient axis
subplot(row, line, estTime - 1);
s = surf(X, Y, Z);
s.FaceAlpha = 1;
s.EdgeColor = 'flat';
s.Marker = 'none';
title('Correlation Coefficient');
xlabel('\it f /\rm Hz');
ylabel('\it \phi /\rm rad');
zlabel('\it R \rm(\itf,\phi\rm)');
set(gca, 'FontSize', 16);
grid on;

% Find maximal value of correlation coefficient
rTemp = R;
[m, idxM] = max(rTemp);
[m1, idxM1] = max(m);
iHat = idxM(idxM1) - 1;
jHat = idxM1 - 1;
rMax = m1;
fHat = fHead + iHat * fInc;
pHat = pHead + jHat * pInc;
fprintf('Estimation No.%d:\n', estTime);
fprintf('Precision: %.5f Hz, %.5f pi rad\n', fInc, pInc/pi);
fprintf('Maximal correlation coefficient: %.5f\n', rMax);
fprintf('Frequency: %f Hz\n', fHat);
fprintf('Phase: %.5f rad (%.5fpi)\n', pHat, pHat/pi);

%% Estimation process 3
% Frequency scale: 0.x(y-1) ~ 0.x(y+1) Hz
% Precision: 0.0002 Hz
% Phase scale: 0.x(y-1)*pi ~ 0.x(y+1)*pi rad
% Phase precision: pi/5000 rad
estTime = estTime + 1;
fe = round(fHat * 100);
pe = round(pHat / pi * 100);

% Sweep freqeuncy and initial phase setting
a0 = 1;                                     % Amplitude (V)
fHead = max(0, (fe - 1)) / 100;             % Starting frequency (Hz)
fEnd = (fe + 1) / 100;                      % Ending frequency (Hz)
fInc = 0.0002;                              % Frequency increment (Hz)
pHead = (pe - 1) / 100 * pi;                % Starting phase (rad)
pEnd = (pe + 1) / 100 * pi;                 % Ending phase (rad)
pInc = pi / 4999;                           % Phase increment (rad)

[R, iMax, jMax] = correlation_sweep_time(a0, fHead, fEnd, fInc, pHead, pEnd, pInc, ...
    idx, xn, Fs);

% Plot correlation coefficient diagram
fn = 0 : iMax;                              % Index matrix of frequency
pn = 0 : jMax;                              % Index matrix of phase
fIdx = fHead + fn * fInc;                   % X coordinate: frequency axis
pIdx = pHead + pn * pInc;                   % Y coordinate: phase axis
[X, Y] = meshgrid(fIdx, pIdx);              % Generate coordinate grid
Z = (R(fn + 1, pn + 1)).';                  % Z coordinate: cross correlation coefficient axis
subplot(row, line, estTime - 1);
s = surf(X, Y, Z);
s.FaceAlpha = 1;
s.EdgeColor = 'flat';
s.Marker = 'none';
title('Correlation Coefficient');
xlabel('\it f /\rm Hz');
ylabel('\it \phi /\rm rad');
zlabel('\it R \rm(\itf,\phi\rm)');
set(gca, 'FontSize', 16);
grid on;

% Find maximal value of correlation coefficient
rTemp = R;
[m, idxM] = max(rTemp);
[m1, idxM1] = max(m);
iHat = idxM(idxM1) - 1;
jHat = idxM1 - 1;
rMax = m1;
fHat = fHead + iHat * fInc;
pHat = pHead + jHat * pInc;
fprintf('Estimation No.%d:\n', estTime);
fprintf('Precision: %.5f Hz, %.5f pi rad\n', fInc, pInc / pi);
fprintf('Maximal correlation coefficient: %.5f\n', rMax);
fprintf('Frequency: %.5f Hz\n', fHat);
fprintf('Phase: %.5f rad (%.5fpi)\n', pHat, pHat / pi);

%% Estimation process 4
% Frequency scale: 0.xy(z-1) ~ 0.xy(z+1) Hz
% Precision: 0.00002 Hz
% Phase scale: 0.xy(z-1)*pi ~ 0.xy(z+1)*pi rad
% Phase precision: pi/500000 rad
estTime = estTime + 1;
fe = round(fHat * 1000);
pe = round(pHat / pi * 1000);

% Sweep freqeuncy and initial phase setting
a0 = 1;                                     % Amplitude (V)
fHead = max(0, (fe - 1)) / 1000;            % Starting frequency (Hz)
fEnd = (fe + 1) / 1000;                     % Ending frequency (Hz)
fInc = 0.00002;                             % Frequency increment (Hz)
pHead = (pe - 1) / 1000 * pi;               % Starting phase (rad)
pEnd = (pe + 1) / 1000 * pi;                % Ending phase (rad)
pInc = pi / 49999;                          % Phase increment (rad)

[R, iMax, jMax] = correlation_sweep_time(a0, fHead, fEnd, fInc, pHead, pEnd, pInc, ...
    idx, xn, Fs);

% Plot correlation coefficient diagram
fn = 0 : iMax;                              % Index matrix of frequency
pn = 0 : jMax;                              % Index matrix of phase
fIdx = fHead + fn * fInc;                   % X coordinate: frequency axis
pIdx = pHead + pn * pInc;                   % Y coordinate: phase axis
[X, Y] = meshgrid(fIdx, pIdx);              % Generate coordinate grid
Z = (R(fn + 1, pn + 1)).';                  % Z coordinate: cross correlation coefficient axis
subplot(row, line, estTime - 1);
s = surf(X, Y, Z);
s.FaceAlpha = 1;
s.EdgeColor = 'flat';
s.Marker = 'none';
title('Correlation Coefficient');
xlabel('\it f /\rm Hz');
ylabel('\it \phi /\rm rad');
zlabel('\it R \rm(\itf,\phi\rm)');
set(gca, 'FontSize', 16);
grid on;

% Find maximal value of correlation coefficient
rTemp = R;
[m, idxM] = max(rTemp);
[m1, idxM1] = max(m);
iHat = idxM(idxM1) - 1;
jHat = idxM1 - 1;
rMax = m1;
fHat = fHead + iHat * fInc;
pHat = pHead + jHat * pInc;
fprintf('Estimation No.%d:\n', estTime);
fprintf('Precision: %.5f Hz, %.5f pi rad\n', fInc, pInc/pi);
fprintf('Maximal correlation coefficient: %.5f\n', rMax);
fprintf('Frequency: %.5f Hz\n', fHat);
fprintf('Phase: %.5f rad (%.5fpi)\n', pHat, pHat/pi);

%% Calculate efficiency indicators
fprintf('\n-------------------------- Display estimation results --------------------------');
totTime = toc;
fPrec = 1 - abs(fHat - ft) / ft;
pPrec = 1 - abs(pHat - pt) / abs(pt);

fprintf("\n");
fprintf('True value: %.5f Hz, %.5f rad (%.5fpi) s\n', ft, pt, pt / pi);
fprintf('Estimation result: %.5f Hz, %.5f rad (%.5fpi) s\n', fHat, pHat, pHat / pi);
fprintf('Sampling time: %.5f s\n', Tt);
fprintf('Estimation time: %.5f s\n', totTime);
fprintf('Total time: %.5f s\n', totTime + Tt);
fprintf('Frequency estimation precision: %.5f%%\n', fPrec * 100);
fprintf('Initial phase estimation precision: %.5f%%\n', pPrec * 100);