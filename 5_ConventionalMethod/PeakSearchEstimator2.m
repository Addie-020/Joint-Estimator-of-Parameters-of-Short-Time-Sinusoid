function [xBest, yBest, info] = PeakSearchEstimator2(xn, Fs)

%
% 
% 
% Input arguments:
%   @xn     : Signal to be estimated
%   @Fs     : Sampling rate
%
% Output arguments:
%   @xBest  : Optimal point (variable)
%   @fBest  : Optimal value of object function
%   @tTot   : Total time of computation
%   @info   : Information of the optimization process
%   @dataLog: Data log of each iteration
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Oct 28, 2022
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

end