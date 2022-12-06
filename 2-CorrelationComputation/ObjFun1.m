function Y = ObjFun1(X, tn, Fs)
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
L = length(tn) - 2;                     % Compute signal length
Xt = (0 : L - 1) / Fs;                  % Time index of samples

Freq = X(1, :);                         % Frequency of current iteration
Phi = X(2, :);                          % Phase of current iteration

% Vecterize settings
F = Freq.';                             % Frequency of current iteration: Rx1
P = repmat(Phi.', 1, L);                % Phase of current iteration: RxL

% Recover test signal and construct estimating signal
Xn = repmat(tn(1 : L), R, 1);           % RxL
Sn = sin(2 * pi * F * Xt + P);          % RxL

% Recover mean and variance of test signal
miuX = repmat(tn(L + 1), R, 1);         % Rx1
sigmaX = repmat(tn(L + 2), R, 1);       % Rx1
% Compute mean and variance of estimating signal
miuS = sum(Sn, 2) / L;                  % Rx1
sigmaS = sqrt(sum((Sn - repmat(miuS, 1, L)).^2, 2) / L);    % Rx1

% Compute cross-correlation coefficient (Person correlation coefficient)
C1 = (Xn - repmat(miuX, 1, L)) ./ repmat(sigmaX, 1, L);     % RxL
C2 = (Sn - repmat(miuS, 1, L)) ./ repmat(sigmaS, 1, L);     % RxL
Rou = sum(C1 .* C2, 2) / (L - 1);       % Rx1

% Compute objective function value
Y = Rou.';

end