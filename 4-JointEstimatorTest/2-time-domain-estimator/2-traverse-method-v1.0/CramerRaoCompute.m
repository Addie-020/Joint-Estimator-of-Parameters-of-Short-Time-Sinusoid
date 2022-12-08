function [varLbAmp, varLbFreq, varLbPha] = CramerRaoCompute(Fs, A, ...
    sigma, N)
%
% Computation of Cramer-Rao Lower Bound (CRLB)
%
% Input arguments:
%   @Fs    : Sampling rate (Hz)
%   @A     : Signal amplitude (A)
%   @sigma : Noise standard variance
%   @N     : Number of sampling points
%
% Output arguments:
%   @varLb : Lower bound of estimation variance
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Oct 25, 2022
%

% Convert SNR to unit form
rou = A^2./(2*sigma.^2);

% Calculate CRLB of amplitude estimation
varLbAmp = 2*sigma.^2/N;

% Calculate CRLB of frequency estimation
freqVar1 = 12*Fs^2;
freqVar2 = (2*pi)^2*rou*N*(N^2-1);
varLbFreq = freqVar1./freqVar2;

% Calculate CRLB of phase estimation
phaVar1 = 2*(2*N-1);
phaVar2 = rou*N*(N+1);
varLbPha = phaVar1./phaVar2;

end