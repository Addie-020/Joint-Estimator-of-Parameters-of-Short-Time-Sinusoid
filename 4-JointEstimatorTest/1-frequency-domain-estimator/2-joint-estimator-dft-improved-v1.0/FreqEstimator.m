function fPeak = FreqEstimator(sigTest, Fs)

%
% Frequency estimator with DFT peak search
% Dichotomous method
% 
% Input arguments:
%   @sigTest: Signal to be estimated
%   @Fs     : Sampling rate
%
% Output arguments:
%   @fPeak: Frequency of the DFT peak
%
% Author        : Zhiyu Shen @Nanjing University
% Establish Date: Dec 5, 2022
% Revised Data  : Dec 5, 2022
%

%%% Step1: Rough Estimation

% Calculate the DFT of the test signal
nSam = length(sigTest);
idxT = 0:nSam-1; 
nFFT = 10*nSam;
testFFT = fft(sigTest, nFFT);
Tw = abs(testFFT);

% Find the peak of DFT
[~, idxPeak] = max(Tw);
p1 = angle(testFFT(idxPeak));

% Compute the frequency value of the peak and the secondary peak
f1 = (idxPeak-2)*Fs/nFFT;
f2 = idxPeak*Fs/nFFT;


%%% Step2: Precise Estimation

% Compute the correlation coefficient at 
% the left end of the frequency search range
b1 = cos(2*pi*f1*idxT+p1);
B1 = abs(fft(b1,nFFT));
r1Vex = corrcoef(Tw,B1);
r1 = min(r1Vex(1,2));

% Compute the correlation coefficient at 
% the right end of the frequency search range
b2 = cos(2*pi*f2*idxT+p1);
B2 = abs(fft(b2,nFFT));
r2Vex = corrcoef(Tw,B2);
r2 = min(r2Vex(1,2));

% Find the maximum correlation coefficient with dichotomous method
while f2-f1 > 0.000000001
    fm = (f1+f2)/2;
    bm = cos(2*pi*fm*idxT+p1);
    Bm = abs(fft(bm,nFFT));
    rmVex = corrcoef(Tw,Bm);
    rm = min(rmVex(1,2));
    if r2 > r1 
        f1 = fm;
        r1 = rm;
    else
        f2 = fm;
        r2 = rm;
    end
end
fPeak = (f1+f2)/2;

end