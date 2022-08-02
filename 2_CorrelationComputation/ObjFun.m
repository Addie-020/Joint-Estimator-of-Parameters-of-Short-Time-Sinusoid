function y = ObjFun(x, tn, fs)
%
% Computation of objective function value
% Objective function is based cross correlation coefficient
%
% Input arguments:
%   @x  : variables (including frequency and phase component)
%   @tn : Sequence to be estimated
%   @fs : Sampling rate
%
% Output arguments:
%   @y  : Objective function value of input variable
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Aug 2, 2022
%

% Set parameters
L = length(tn);                 % Compute signal length
xt = (0 : L - 1) / fs;          % Time index of samples

f = x(1);                       % Frequency of current iteration
phi = x(2);                     % Phase of current iteration

% Construct estimating signal
sn = sin(2 * pi * f * xt + phi);

% Compute objective function value
c = corrcorf(tn, sn);
y = exp(5 * c(2));

end