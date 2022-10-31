function varLb = CramerRaoCompute(Fs, snrSig, N)
%
% Computation of Cramer-Rao Lower Bound (CRLB)
%
% Input arguments:
%   @Fs    : Sampling rate
%   @snrSig: SNR
%   @N     : Number of sampling points
%
% Output arguments:
%   @varLb : Lower bound of estimation variance
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Oct 25, 2022
%

rou = 10.^(snrSig/10);
tempVar1 = 12 * Fs^2;
tempVar2 = (2*pi)^2 * rou * N^3;

varLb = tempVar1 ./ tempVar2;

end