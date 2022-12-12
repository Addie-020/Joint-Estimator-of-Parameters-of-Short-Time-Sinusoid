% Description:  Plot MSE-SNR curve of Joint Estimator
% Projet:       Joint Estimatior of Frequency and Phase
% Date:         Dec 10, 2022
% Author:       Zhiyu Shen

close all
clear
clc

%% Fetch Data

% Estimation condition
Fs = 5;
Ns = 32;

% Read Data from Excel Table
dataIn = readmatrix(['E:\1-academic\2-projects\1-2-short-signal' ...
    '-estimation\3-working\01-test-1211\time-MSE-SNR-options.xlsx'], ...
    'Sheet', 'pso', 'Range', 'A4');

% Assign data to corresponding vector
SNRdB = dataIn(:, 1);
mseLbFreq = dataIn(:, 2);
mseLbPhas = dataIn(:, 3);
% MSE vector: store MSE under diffrent tolFunVal of PSO in column
% 1: 1e-15    2: 1e-12    3: 1e-9    4: 1e-6
mseFreq = dataIn(:, [4 6 8 10]);
msePhas = dataIn(:, [5 7 9 11]);


%% Plot

% Plot frequency MSE-SNR curve
fErrPlt = figure(1);
fErrPlt.Name = "Relationship between frequency MSE and SNR";
fErrPlt.WindowState = 'maximized';
semilogy(SNRdB, mseLbFreq, 'LineWidth', 2, 'Color', '#77AC30', ...
    'Marker', 'square', 'LineStyle', '-.');
hold on
semilogy(SNRdB, mseFreq(:,1), 'LineWidth', 2, 'Color', '#A2142F', ...
    'Marker', 'x', 'LineStyle', ':');
semilogy(SNRdB, mseFreq(:,2), 'LineWidth', 2, 'Color', '#7E2F8E', ...
    'Marker', '*', 'LineStyle', ':');
semilogy(SNRdB, mseFreq(:,3), 'LineWidth', 2, 'Color', '#EDB120', ...
    'Marker', 'o', 'LineStyle', ':');
semilogy(SNRdB, mseFreq(:,4), 'LineWidth', 2, 'Color', '#0072BD', ...
    'Marker', '.', 'LineStyle', ':');
hold off
xlabel("SNR (dB)", "Interpreter", "latex");
ylabel("$MSE_{frequency}$", "Interpreter", "latex");
legend("CRLB", "1e-15", "1e-12", "1e-9", "1e-6");
set(gca, 'Fontsize', 20);

% Plot phase MSE-SNR curve
pErrPlt = figure(2);
pErrPlt.Name = "Relationship between frequency MSE and SNR";
pErrPlt.WindowState = 'maximized';
semilogy(SNRdB, mseLbPhas, 'LineWidth', 2, 'Color', '#77AC30', ...
    'Marker', 'square', 'LineStyle', '-.');
hold on
semilogy(SNRdB, msePhas(:,1), 'LineWidth', 2, 'Color', '#A2142F', ...
    'Marker', 'x', 'LineStyle', ':');
semilogy(SNRdB, msePhas(:,2), 'LineWidth', 2, 'Color', '#7E2F8E', ...
    'Marker', '*', 'LineStyle', ':');
semilogy(SNRdB, msePhas(:,3), 'LineWidth', 2, 'Color', '#EDB120', ...
    'Marker', 'o', 'LineStyle', ':');
semilogy(SNRdB, msePhas(:,4), 'LineWidth', 2, 'Color', '#0072BD', ...
    'Marker', '.', 'LineStyle', ':');
hold off
xlabel("SNR (dB)", "Interpreter", "latex");
ylabel("$MSE_{phase}$", "Interpreter", "latex");
legend("CRLB", "1e-15", "1e-12", "1e-9", "1e-6");
set(gca, 'Fontsize', 20);


%% Fit the Curve with Linear Regression

% Extract frequency MSE data from matrix
mseFreq15 = log10(mseFreq(:,1));
mseFreq12 = log10(mseFreq(:,2));
mseFreq09 = log10(mseFreq(:,3));
mseFreq06 = log10(mseFreq(:,4));

% Extract phase MSE data from matrix
msePhas15 = log10(msePhas(:,1));
msePhas12 = log10(msePhas(:,2));
msePhas09 = log10(msePhas(:,3));
msePhas06 = log10(msePhas(:,4));

[mm, nn] = size(mseFreq);
mseFreqFit = zeros(mm, nn);
freqParam15 = polyfit(SNRdB, mseFreq15, 2);
mseFreqFit(:,1) = polyval(freqParam15, SNRdB);
