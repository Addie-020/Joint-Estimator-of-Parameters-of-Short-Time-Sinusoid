function [freqMse, phaMse, timeMean, timeVar] = PeakSearchTest(xn, ...
    ft, pt, Fs, Tt, numEst)

%
% Test function of Joint Estimator
% Run estimator for multiple times
% 
% Input arguments:
%   @xn     : Signal to be estimated
%   @ft     : Frequency of signal to be estimated
%   @pt     : Phase of signal to be estimated
%   @Fs     : Sampling rate (Hz)
%   @Tt     : Total time of sampling (s)
%   @numEst : Estimation times for each test
%   @maxIter: Maximum iteration time for each estimation
%
% Output arguments:
%   @freqMse : MSE of frequency estimated
%   @phaMse  : MSE of phase estimated
%   @timeMean: Mean of time of each estimation
%   @timeVar : Variance of time of each estimation
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Nov 1, 2022
%

%%% Estimation Process

% Estimate loop
timeTot = zeros(1, numEst);         % Estimation time for each iteration
fe = zeros(1, numEst);              % Estimated frequency of each iteration
pe = zeros(1, numEst);              % Estimated phase of each iteration
for i = 1 : numEst
    tic
    [xBest, ~] = PeakSearchEstimator2(xn, Fs);
    timeTot(i) = toc;
    % Assign results
    fe(i) = xBest(1);
    pe(i) = xBest(2);
end


%%% Process result

% Process estimation time
timeEst = timeTot + Tt;
timeMean = sum(timeEst) ./ numEst;
timeVar = sum((timeEst-timeMean).^2) / numEst;

% Calculate error
freqMse = sum((fe-ft).^2) ./ numEst;
phaMse = sum((pe-pt).^2) ./ numEst;

end % end: function PeakSearchTest



%%%% Function "PeakSearchEstimator2"

function [xBest, yBest] = PeakSearchEstimator2(xn, Fs)

%
% Estimate frequency, amplitude and phase
% Search the peak of DTFT spectrum (Improved with CZT)
% 
% Input arguments:
%   @xn     : Signal to be estimated
%   @Fs     : Sampling rate (Hz)
%
% Output arguments:
%   @xBest  : Optimal point (variable)
%   @fBest  : Optimal value of object function
%


%%% Pre-Processing
Ns = length(xn);                            % Number of samples
xn = xn - sum(xn)/Ns;                       % De-mean


% Calculate amplitude spectrum
xw = fft(xn)*2/Ns;
xwAmp = abs(xw);

% Plot amplitude spectrum
subplot(2, 1, 1);
plot((0:Ns/2-1)*Fs/Ns, xwAmp(1:Ns/2));
title("Signal's Amplitude Spectrum", "Interpreter", "latex");
xlabel("Frequency $f\ (Hz)$", "Interpreter", "latex");
ylabel("Amplitude $A\ (V)$", "Interpreter", "latex");
set(gca, 'Fontsize', 20);


%%% Find the Peak of DFT

% Find peak index
% If the index is not unique, then use the first one
idxPeak = find(xwAmp==max(xwAmp(1:Ns/2)));
idxPeak = idxPeak(1);

% Calculate peak frequency according to index
f0 = (idxPeak-1)*Fs/Ns;


%%% CZT

% Define parameters
M = 10000;                      % Variable controls the frequency resolution
fLb = f0 - Fs/Ns;               % Lowe bound of search frequency
fInc = 2*Fs/(Ns*M);             % Search frequency resolution
fSearch = fLb : fInc : fLb+(M-1)*fInc;

% CZT calculation
xSum = xn(1).*ones(1, length(fSearch));
for n = 2 : Ns
    xSum = xSum + xn(n)*exp(-1i*2*pi*fSearch*(n-1)/Fs);
end
Y = abs(xSum)*2/Ns;

% Plot CZT result
subplot(2, 1, 2);
plot(fSearch, Y);
title("Signal's CZT Result", "Interpreter", "latex");
xlabel("Frequency $f\ (Hz)$", "Interpreter", "latex");
ylabel("Amplitude $A\ (V)$", "Interpreter", "latex");
set(gca, 'Fontsize', 20);


%%% Peak Search

% Find peak index
idxPeak = find(Y==max(Y));
idxPeak = idxPeak(1);

% Calculate peak frequency
fPeak = fLb + (idxPeak-1)*fInc;
YPeak = Y(idxPeak);

% Calculate peak phase
xwPha = angle(xSum);
xwPha = unwrap(xwPha);
pPeak = xwPha(idxPeak);

% Correct phase
while 0 < 1
    if pPeak > 2*pi
        pPeak = pPeak - 2*pi;
    elseif pPeak < 0
        pPeak = pPeak + 2*pi;
    else
        break;
    end
end


%%% Output
xBest = [fPeak, pPeak];
yBest = YPeak;

end % end: function PeakSearchEstimator2