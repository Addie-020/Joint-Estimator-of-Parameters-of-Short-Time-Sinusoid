function [xn, L] = WaveGen(at, ft, pt, Fs, Tt)
%
% Generate a test sinusoid sequence according to input parameters
% 
% Input arguments:
%   @at: Amplitude of sinusoid
%   @ft: Frequency of sinusoid
%   @pt: Initial phase of sinusoid
%   @Fs: Sampling frequency
%   @Tt: Total sample time
%
% Output arguments:
%   @xn : Output sinusoid sequence
%   @L  : Length of output sequence
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Aug 3, 2022
%

L = Tt*Fs;                          % Total sampling points
tIdx = (0:L-1)/Fs;                  % Time index
xn = at*cos(2*pi*ft*tIdx+pt);       % Test signal

end