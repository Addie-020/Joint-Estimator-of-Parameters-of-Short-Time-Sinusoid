function Y = ObjFun(X, Ct, Fs)
%
% Computation of objective function value
% Objective function is based cross correlation coefficient
% X is a nParticles*nvars vector, dimension is shown in row
%
% Input arguments:
%   @X  : variables (including frequency and phase component)
%   @Ct : Covariance information of sequence to be estimated
%   @Fs : Sampling rate
%
% Output arguments:
%   @y  : Objective function value of input variable
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Sept 15, 2022
%

% Input vector size validation
[~, m] = size(X);
if m ~= 2
    error('X is not of a valid size!')
end % end: if

% Set parameters
nSequence = length(Ct);                             % Compute signal length
xt = (0 : nSequence-1) / Fs;                        % Time index of samples

% Set freuqency and phase vector
freq = X(:, 1);                                     % nParticles*1
phi = X(:, 2);                                      % nParticles*1

% Construct estimating signal
Sn = sin(2*pi*freq*xt + phi);                       % nParticles*nSequence

% Compute mean and variance of estimating signal
miuS = sum(Sn, 2) / nSequence;                      % nParticles*1
sigmaS = sqrt(sum((Sn - miuS).^2, 2) / nSequence);  % nParticles*1

% Compute cross-correlation coefficient (Person correlation coefficient)
Ce = (Sn - miuS) ./ sigmaS;                         % nParticles*nSequence
Rou = Ce * Ct.' / (nSequence - 1);                  % nParticles*1

% Compute objective function value
Y = 8 - exp(Rou+1);                                 % nParticles*1

end