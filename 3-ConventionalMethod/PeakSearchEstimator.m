function [xBest, yBest] = PeakSearchEstimator(xn, Fs)

%
% Estimate frequency, amplitude and phase
% Search the peak of DTFT spectrum
% 
% Input arguments:
%   @xn     : Signal to be estimated
%   @Fs     : Sampling rate
%
% Output arguments:
%   @xBest  : Optimal point (variable)
%   @fBest  : Optimal value of object function
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Oct 28, 2022
%


%%% Pre-Processing
Ns = length(xn);                            % Number of samples
xn = xn - sum(xn)/Ns;                       % De-mean


%%% DFT

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

% Calculate phase spectrum
xwPha = angle(xw);
xwPha = unwrap(xwPha);

% Plot phase spectrum
subplot(2, 1, 2);
plot((0:Ns/2-1)*Fs/Ns, xwPha(1:Ns/2));
title("Signal's Phase Spectrum", "Interpreter", "latex");
xlabel("Frequency $f\ (Hz)$", "Interpreter", "latex");
ylabel("Phase $\phi (rad)$", "Interpreter", "latex");
set(gca, 'Fontsize', 20);


%%% Find the Peak of DFT

% Find peak index
% If the index is not unique, then use the first one
idxPeak = find(xwAmp==max(xwAmp(1:Ns/2)));
idxPeak = idxPeak(1);

% Calculate peak frequency, phase and amplitude according to index
fPeak = (idxPeak-1)*Fs/Ns;
aPeak = xwAmp(idxPeak);
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
yBest = aPeak;

end