function Y = ObjFunComp(X, tn, Fs)
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
R = size(X, 2);                         % Repetition of input variables
L = length(tn);                         % Compute signal length
Xt = (0 : L - 1) / Fs;                  % Time index of samples

F = X(1, :);                            % Frequency of current iteration
P = X(2, :);                            % Phase of current iteration

% Recover test signal and construct estimating signal
Xn = tn(1 : L);

Rou = zeros(1, R);
for i = 1 : R
    Sn = sin(2 * pi * F(i) * Xt + P(i));
    C = corrcoef(Xn, Sn);
    Rou(i) = C(2);
end

% Compute objective function value
Y = Rou;

end