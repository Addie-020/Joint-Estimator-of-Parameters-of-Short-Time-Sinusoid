function xBest = MatchedSpectrum(xn, Fs)
%
% Frequency estimator with DFT peak search
% Dichotomous method
% 
% Input arguments:
%   @xn: Signal to be estimated
%   @Fs: Sampling frequency
%
% Output arguments:
%   @xPeak: Frequency and phase result of the DFT peak search
%
% Author         : Zhiyu Shen @Nanjing University
% Estabilish date: Dec 8, 2022
% Revise data    : Dec 9, 2022
%

% Add window and perform DFT on the test signal
nSam = length(xn);                              % Number of original signal samples
nFFT = 2^nextpow2(nSam);                        % Number of FFT points
winSig = ones(1,nFFT);                          % Window signal (Rectangular Window)
% idxWin = 0 : 1 : nFFT-1;
% winSig = 0.54 - 0.46*cos(2*pi*idxWin/nFFT);     % Window signal (Hamming Window)
xnWin = [xn, zeros(1,nFFT-nSam)].*winSig;       % Zero padding and add window
xnFFT = fft(xnWin,nFFT);

% Compute the frequency spectrum of test signal
Xn1 = abs(xnFFT/nFFT);
Xw = Xn1(1:nFFT/2);
Xw(2:end-1) = 2*Xw(2:end-1);

% Compute mean and variance of test signal
miu0 = sum(Xw)/nFFT;
sigma0 = sqrt(sum((Xw-miu0).^2)/nFFT);

% Compute signal information for correlation computation
Ct = (Xw-miu0) ./ sigma0;
Ct = [Ct, nSam, nFFT];

% Find the peak of DFT
[~, idxPeak] = max(Xw);
p1 = angle(xnFFT(idxPeak));

% Compute the frequency value of the peak and the secondary peak
f1 = (idxPeak-2)*Fs/nFFT;
f2 = idxPeak*Fs/nFFT;

% Compute the correlation coefficient at 
% the left and right end of the frequency search range
r1 = ObjFunFreq([f1 p1], Ct, Fs); 
r2 = ObjFunFreq([f2 p1], Ct, Fs); 

% Find the maximum correlation coefficient with dichotomous method
while f2-f1 > 1e-9
    fm = (f1+f2)/2;
    rm = ObjFunFreq([fm p1], Ct, Fs); 
    if r2 < r1 
        f1 = fm;
        r1 = rm;
    else
        f2 = fm;
        r2 = rm;
    end
end
fPeak = (f1+f2)/2;

% Estimate phase roughly
if idxPeak == 1
    f1 = (idxPeak-1)*Fs/nFFT;
    f2 = idxPeak*Fs/nFFT;
    idx1 = idxPeak;
    idx2 = idxPeak+1;
else
    if Xw(idxPeak-1) > Xw(idxPeak+1)
        f1 = (idxPeak-2)*Fs/nFFT;
        f2 = (idxPeak-1)*Fs/nFFT;
        idx1 = idxPeak-1;
        idx2 = idxPeak;
    else
        f1 = (idxPeak-1)*Fs/nFFT;
        f2 = idxPeak*Fs/nFFT;
        idx1 = idxPeak;
        idx2 = idxPeak+1;
    end
end

% Compute the angle value at peak and correct it
p1 = angle(xnFFT(idx1));
p2 = angle(xnFFT(idx2));
if abs(p1-p2) > 5
    if p1 > p2
        p2 = p2 + 2*pi;
        p0 = angle(xnFFT(idx2)) + (fPeak-f2)*(p1-p2)/(f1-f2);
    else
        p1 = p1 + 2*pi;
        p0 = angle(xnFFT(idx2)) + (fPeak-f2)*(p1-p2)/(f1-f2);
    end
else
    p0 = p2 + (fPeak-f2)*(p1-p2)/(f1-f2);
end
pPeak = mod(p0, 2*pi);

% Generate output vector
xBest = [fPeak, pPeak];

end % end: function RoughEstimate



%%%% Function "ObjFunFreq"

function Y = ObjFunFreq(var, Ct, Fs)
%
% Computation of objective function value
% Objective function is based cross correlation coefficient
% X is a nParticles*nvars vector, dimension is shown in row
%
% Input arguments:
%   @var: variables (including frequency and phase component)
%   @Ct : Covariance information of test signal's frequency spectrum
%   @Fs : Sampling rate
%
% Output arguments:
%   @y  : Objective function value of input variable
%

% Input vector size validation
[nSig, nVar] = size(var);
if nVar ~= 2
    error('X is not of a valid size!')
end % end: if

% Set parameters
nFFT = Ct(end);                                 % Fetch FFT length
nSam = Ct(end-1);                               % Fetch original signal length
Ct = Ct(1:end-2);                               % Fetch real values of Ct
idxT = (0:nSam-1)/Fs;                           % Time index of samples

% Set freuqency and phase vector
varF = var(:,1);                                % nParticles*1
varP = var(:,2);                                % nParticles*1

% Construct estimating signal
sn = cos(2*pi*varF*idxT+varP);                  % nParticles*NFFT

% Add window and perform DFT on the constructed signal
winSig = ones(1,nFFT);                          % Window signal (Rectangular Window)
% idxWin = 0:1:nFFT-1;
% winSig = 0.54-0.46*cos(2*pi*idxWin/nFFT);       % Window signal (Hamming Window)
snWin = [sn, zeros(nSig,nFFT-nSam)].* ...
    repmat(winSig,nSig,1);                      % Zero padding and add window
snFFT = fft(snWin,nFFT,2);

% Compute the frequency spectrum of constructed signal
Sn1 = abs(snFFT./nFFT);
Sn = Sn1(:,1:nFFT/2);
Sn(:,2:end-1) = 2*Sn(:,2:end-1);

% Compute mean and variance of estimating signal
miuS = sum(Sn,2)/nFFT;                          % nParticles*1
sigmaS = sqrt(sum((Sn-miuS).^2,2)/nFFT);        % nParticles*1

% Compute cross-correlation coefficient (Person correlation coefficient)
Ce = (Sn-miuS)./sigmaS;                         % nParticles*nSequence
Rou = Ce*Ct.'/(nFFT-1);                         % nParticles*1

% Compute objective function value
Y = 8-exp(Rou+1);                               % nParticles*1

end % end: function ObjFunFreq