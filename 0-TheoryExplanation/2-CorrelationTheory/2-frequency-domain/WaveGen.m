function [xn, Ns] = WaveGen(at, ft, pt, Fs, Ns)
%
% Generate a test sinusoid sequence according to input parameters
% 
% Input arguments:
%   @at: Amplitude of sinusoid
%   @ft: Frequency of sinusoid
%   @pt: Initial phase of sinusoid
%   @Fs: Sampling frequency
%   @Ns: Number of samples
%
% Output arguments:
%   @xn : Output sinusoid sequence
%   @L  : Length of output sequence
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Aug 3, 2022
%

tIdx = (0:Ns-1)/Fs;                  % Time index
xn = at*cos(2*pi*ft*tIdx+pt);       % Test signal

end