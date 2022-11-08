function [freqMse, phaMse, timeMean, timeVar] = PhaseDiffTest(xn, ...
    ft, pt, Fs, Tt, numEst)

%
% Test function of Phase Difference Method
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
    [xBest, ~] = PhaseDiff(xn, Fs);
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

end % end: function PhaseDiffTest



%%%% Function "PhaseDiff"

function [xBest, yBest] = PhaseDiff(xn, Fs)

%
% Estimate frequency, amplitude and phase
% Search the peak of DTFT spectrum (Improved with phase difference)
%
% Input arguments:
%   @xn: Signal to be estimated
%   @Fs: Sampling rate
%
% Output arguments:
%   @xBest: Optimal point (variable)
%   @fBest: Optimal value of object function
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Nov 3, 2022
%

%%% Pre-Processing

% Define parameters
Ns = length(xn);            % Number of samples
Ts = 1/Fs;                  % Sampling interval (s)
N2 = 2*round(Ns/4);
N1 = Ns - N2;

% Process sequence
xn = xn - sum(xn)/Ns;       % De-mean
win = hanning(Ns).';        % Generate a window function
xn = xn.*win;               % Add window to signal
x1 = xn;
x2 = xn(N1+1 : end);


%%% Find the Peak with DFT

% Calculate amplitude spectrum
xw = fft(xn)*2/Ns;
xwAmp = abs(xw);

% Find peak index
% If the index is not unique, then use the first one
idxPeak = find(xwAmp==max(xwAmp(1:Ns/2)));
idxPeak = idxPeak(1);

% Calculate peak frequency according to index
f0 = (idxPeak-1)*Fs/Ns;
f0Lb = f0 - 0.1*Fs/Ns;          % Lowe bound of search frequency


%%% Improve the Accuracy with Phase Difference

% Calculate two phases
f = f0Lb;
y1 = DTFT_f(xn, f, Ns, Fs);
y2 = DTFT_f(xn, f, N2, Fs);
xwPha = phase(y1);
xwPha = unwrap(xwPha);
pTemp = xwPha;
while 0 < 1 
    if pTemp > 2*pi
        pTemp = pTemp - 2*pi;
    elseif pTemp<0
        pTemp = pTemp + 2*pi;
    else
        break
    end
end
pha1=pTemp;
xwPha = phase(y2);
xwPha = unwrap(xwPha);
pTemp = xwPha;
while 0 < 1
    if pTemp > 2*pi
        pTemp = pTemp - 2*pi;
    elseif pTemp < 0
        pTemp = pTemp + 2*pi;
    else
        break
    end
end
pha2 = pTemp;

% Correct the result
fPeak = f - (pha2-pha1)*2/(N2*pi*Ts);
y = DTFT_f(xn, fPeak, Ns, Fs);
aPeak = abs(y);
xwPha = phase(y);
xwPha = unwrap(xwPha);
pTemp = xwPha;
while 0 < 1
    if pTemp > 2*pi
        pTemp = pTemp - 2*pi;
    elseif pTemp < 0
        pTemp = pTemp + 2*pi;
    else
        break
    end
end
pPeak=pTemp;


%%% Output
xBest = [fPeak, pPeak];
yBest = aPeak;

end % end: function PhaseDiff



%%%% Function "DTFT_f"

function y = DTFT_f(xn, f0, Ns, Fs)

%
% Calculate DTFT manually for single point
%
% Input arguments:
%   @xn: Input sequence to be processed
%   @f0: Frequency of the point
%   @Ns: Number of sampled ponits
%   @Fs: Sampling rate (Hz)
%
% Output arguments:
%   @y: DTFT value of the input point
%

y = 0;
for n = 0 : Ns-1 
   y = y + xn(n+1)*exp(-1i*2*pi*f0/Fs*n);
end
y = y*4/Ns;

end % end: function DTFT_f