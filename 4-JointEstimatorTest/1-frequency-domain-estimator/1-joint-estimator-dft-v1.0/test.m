close all
clear
clc

% Generate original signal sequence
at = 1;
ft = 0.5;
pt = pi/3;
Fs = 5;
Ns = 32;
xt = (0:Ns-1)/Fs;
xn0 = at*cos(2*pi*ft*xt+pt);
xn = xn0;

% Add window and perform DFT on the test signal
NFFT = 2^nextpow2(Ns);                          % Number of FFT points
idxWin = 0:1:NFFT-1;
winSig = 0.54-0.46*cos(2*pi*idxWin/NFFT);       % Window signal (Hamming Window)
xnWin = [xn, zeros(1,NFFT-Ns)].*winSig;         % Zero padding and add window
xnFFT = fft(xnWin,NFFT);

% Compute the frequency spectrum of test signal
Xn1 = abs(xnFFT/NFFT);
Xn = Xn1(1:NFFT/2);
Xn(2:end-1) = 2*Xn(2:end-1);

% Compute mean and variance of test signal
miu0 = sum(Xn) / NFFT;
sigma0 = sqrt(sum((Xn-miu0).^2) / NFFT);

% Compute signal information for correlation computation
Ct = (Xn-miu0) ./ sigma0;
Ct = [Ct, Ns, NFFT];

% Compute objective function value related to correlation coefficient
var = [0.1 pi/2; 0.5 pi/3; 0.7 pi/5];
Y = ObjFunFreq(var, Ct, Fs);