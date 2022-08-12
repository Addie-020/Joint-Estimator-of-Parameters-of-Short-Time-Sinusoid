function [xn, L, idx] = WaveGen(at, ft, pt, Fs, Tt)

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
%   @idx: Index of output sequence
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Aug 3, 2022
%

L = Tt * Fs;                            % Total sampling points
idx = 0 : L - 1;                        % Index of sampling points
xt = idx / Fs;                          % Time index
xn = at * sin(2 * pi * ft * xt + pt);   % Test signal

end