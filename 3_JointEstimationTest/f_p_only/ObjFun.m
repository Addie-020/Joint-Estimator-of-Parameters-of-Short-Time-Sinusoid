function Y = ObjFun(X, Ct, Fs)
%
% Computation of objective function value
% Objective function is based cross correlation coefficient
%
% Input arguments:
%   @X  : variables (including frequency and phase component)
%   @Ct : Necessary information of sequence to be estimated
%   @Fs : Sampling rate
%
% Output arguments:
%   @y  : Objective function value of input variable
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Aug 2, 2022
%

% Set parameters
R = size(X, 2);                         % Repetition of input variables
L = length(Ct);                         % Compute signal length
Xt = (0 : L - 1) / Fs;                  % Time index of samples

Freq = X(1, :);                         % Frequency of current iteration
Phi = X(2, :);                          % Phase of current iteration

% Vecterize settings
F = Freq.';                             % Frequency of current iteration: Rx1
P = repmat(Phi.', 1, L);                % Phase of current iteration: RxL

% Construct estimating signal
Sn = sin(2 * pi * F * Xt + P);                              % RxL

% Compute mean and variance of estimating signal
miuS = sum(Sn, 2) / L;                                      % Rx1
sigmaS = sqrt(sum((Sn - repmat(miuS, 1, L)).^2, 2) / L);    % Rx1

% Compute cross-correlation coefficient (Person correlation coefficient)
C1 = repmat(Ct, R, 1);                                      % RxL
C2 = (Sn - repmat(miuS, 1, L)) ./ repmat(sigmaS, 1, L);     % RxL
Rou = sum(C1 .* C2, 2) / (L - 1);                           % Rx1

% Compute objective function value
Y = 2 - Rou.';

end