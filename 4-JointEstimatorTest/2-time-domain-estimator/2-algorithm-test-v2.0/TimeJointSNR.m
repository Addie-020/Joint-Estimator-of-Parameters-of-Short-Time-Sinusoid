% Description:  Program for Method Comparison in MSE-SNR performance
% Projet:       Joint Estimatior of Frequency and Phase
% Date:         Feb 5, 2022
% Author:       Zhiyu Shen

clear
close all
clc


%% Set Up Estimation Options

% Set frequency and phase range
fLb = 49;
fUb = 51;
pLb = 0;
pUb = 2*pi;
paramRange = [fLb, fUb, pLb, pUb];

% Set sampling parameters (Hz)
Fs = 50;
% Fs = 10*(fUb+fLb)/2;

% Set sampling points
Ns = 32;                             % Total sampling points
Tt = Ns/Fs;                         % Total time of sampling (s)

% Set noise figure
at = 1;                             % Signal amplitude
SNRdB = 0:5:80;                     % SNR (dB)

% Generate signal time index
xt = (0:Ns-1)/Fs;                   % Time index


%% Estimation Process

% Estimator settings
options.maxIter = 3;                        % Search times for each estimation
optionsPso.tolFunValue = 1e-9;              % Termination tolerance on function value [1e-9]
optionsPso.swarmSize = 100;                 % Number of particles [100, 10nvars]
optionsGrad.Algorithm = 'interior-point';
optionsGrad.ConstraintTolerance = 1e-6;     % [1e-6]
optionsGrad.StepTolerance = 1e-10;          % [1e-10]

% Define estimator options and allocate vector memories
numEst = 1000;                      % Number of estimations
numSNR = length(SNRdB);             % Number of SNR points
numMet = 2;                         % Number of methods
mseFreq = zeros(numSNR, numMet);    % MSE of frequency estimation
msePhas = zeros(numSNR, numMet);    % MSE of phase estimation
rmseFreq = zeros(numSNR, numMet);   % RMSE of frequency estimation
rmsePhas = zeros(numSNR, numMet);   % RMSE of phase estimation
mseLbFreq = zeros(numSNR, 1);       % MSE lower bound of frequency estimation
mseLbPhas = zeros(numSNR, 1);       % MSE lower bound of phase estimation
rmseLbFreq = zeros(numSNR, 1);      % RMSE lower bound of frequency estimation
rmseLbPhas = zeros(numSNR, 1);      % RMSE lower bound of phase estimation
    
% Estimate loop
poolobj = parpool(12);
parfor ii = 1 : numSNR
% for ii = 1 : numSNR
    
    % Allocate memeory space for recording vectors
    errFreq = zeros(numEst, numMet);        % Frequency estimation error vector
    errPhas = zeros(numEst, numMet);        % Phase estimation error vector
    
    % Generate noise sequence
    sigmaN = at / 10.^(SNRdB(ii)/20);       % Standard variance of noise
    sigNoise = sigmaN * randn(1, Ns);       % Additive white Gaussian noise  
    
    % Estimation of single SNR
    for jj = 1 : numEst
        
        % Generate signal sequence and add noise
        ft = fLb + 0.01*randi([0 round(100*(fUb-fLb))]);
        pt = pLb + 0.01*randi([0 round(100*(pUb-pLb))]);
        x0 = at * cos(2*pi*ft*xt + pt);
        xn = x0 + sigNoise;

        % ---------- Joint estimator ----------
        [xBest, ~, ~] = JointEstimatorTime(xn, Fs, paramRange, options, [], []);
        fe = xBest(1);
        pe = xBest(2);
        errFreq(jj,1) = abs(fe-ft);
        errPhas(jj,1) = min(abs([pe-pt; pe-pt+2*pi; pe-pt-2*pi]));

        % ---------- Mao's joint method ----------
        [xBest, ~] = MaoJoint(xn, Fs)
        fe = xBest(1);
        pe = xBest(2);
        errFreq(jj,2) = abs(fe-ft);
        errPhas(jj,2) = min(abs([pe-pt; pe-pt+2*pi; pe-pt-2*pi]));
        
    end % end for

    % Compute RMSE
    mseFreq(ii,:) = sum(errFreq.^2)./numEst;
    msePhas(ii,:) = sum(errPhas.^2)./numEst;
    rmseFreq(ii,:) = sqrt(mseFreq(ii,:));
    rmsePhas(ii,:) = sqrt(msePhas(ii,:));

    % Compute CRLB
    [~, mseLbFreq(ii), mseLbPhas(ii)] = CramerRaoCompute(Fs, at, sigmaN, Ns);
    rmseLbFreq(ii) = sqrt(mseLbFreq(ii));
    rmseLbPhas(ii) = sqrt(mseLbPhas(ii));

    % Print iteration info
    fprintf('Iteration No.%d\n', ii);

end
delete(poolobj);


%% Output

fprintf('\n-------- RMSE-SNR Test of Joint Estimator --------\n');
fprintf('Fs = %.2f Hz\n', Fs);
fprintf('Sampling points = %d\n', Ns);
fprintf('Sampling time = %.2f s\n', Tt);
fprintf('Number of estimations per SNR = %d\n', numEst);


%% Plot Figures

% Plot relationship between MSE and SNR
% errPlt = figure(1);
% errPlt.Name = "Relationship between MSE and SNR";
% errPlt.WindowState = 'maximized';

% Plot frequency MSE-SNR curve
fErrPlt = figure(1);
fErrPlt.Name = "Relationship between frequency MSE and SNR";
fErrPlt.WindowState = 'maximized';
hold on
plot(SNRdB, log10(mseLbFreq), 'LineWidth', 2, 'Color', '#77AC30', ...
    'Marker', 'square', 'LineStyle', '-.');
plot(SNRdB, log10(mseFreq(:,1)), 'LineWidth', 2, 'Color', '#A2142F', ...
    'Marker', 'x', 'LineStyle', '--');
plot(SNRdB, log10(mseFreq(:,2)), 'LineWidth', 2, 'Color', '#7E2F8E', ...
    'Marker', '*', 'LineStyle', '--');
% plot(SNRdB, log10(mseFreq(:,3)), 'LineWidth', 2, 'Color', '#EDB120', ...
%     'Marker', 'o', 'LineStyle', ':');
% plot(SNRdB, log10(mseFreq(:,4)), 'LineWidth', 2, 'Color', '#0072BD', ...
%     'Marker', '.', 'LineStyle', ':');
% plot(SNRdB, log10(mseFreq(:,5)), 'LineWidth', 2, 'Color', '#D95319', ...
%     'Marker', 'diamond', 'LineStyle', ':');
hold off
xlabel("SNR (dB)", "Interpreter", "latex");
ylabel("$\log_{10}(MSE_{frequency})$", "Interpreter", "latex");
legend("CRLB", "Joint(Time)", "Joint(Mao)");
set(gca, 'Fontsize', 20);

% Plot phase MSE-SNR curve
pErrPlt = figure(2);
pErrPlt.Name = "Relationship between frequency MSE and SNR";
pErrPlt.WindowState = 'maximized';
hold on
plot(SNRdB, log10(mseLbPhas), 'LineWidth', 2, 'Color', '#77AC30', ...
    'Marker', 'square', 'LineStyle', '-.');
plot(SNRdB, log10(msePhas(:,1)), 'LineWidth', 2, 'Color', '#A2142F', ...
    'Marker', 'x', 'LineStyle', '--');
plot(SNRdB, log10(msePhas(:,2)), 'LineWidth', 2, 'Color', '#7E2F8E', ...
    'Marker', '*', 'LineStyle', '--');
% plot(SNRdB, log10(msePhas(:,3)), 'LineWidth', 2, 'Color', '#EDB120', ...
%     'Marker', 'o', 'LineStyle', ':');
% plot(SNRdB, log10(msePhas(:,4)), 'LineWidth', 2, 'Color', '#0072BD', ...
%     'Marker', '.', 'LineStyle', ':');
% plot(SNRdB, log10(msePhas(:,5)), 'LineWidth', 2, 'Color', '#D95319', ...
%     'Marker', 'diamond', 'LineStyle', ':');
hold off
xlabel("SNR (dB)", "Interpreter", "latex");
ylabel("$\log_{10}(MSE_{phase})$", "Interpreter", "latex");
legend("CRLB", "Joint(Time)", "Joint(Mao)");
set(gca, 'Fontsize', 20);

