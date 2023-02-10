function Y = CorrCal(Xw, Yw)
%
% Computation of objective function value
% Objective function is based cross correlation coefficient
%
% Input arguments:
%   @Xw : Signal 1
%   @Yw : Signal 2
%   @Fs : Sampling rate
%
% Output arguments:
%   @y  : Objective function value of input variable
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Feb 2, 2023
%

% Set parameters
N = length(Xw);                         % Compute signal length

% Compute statistical characteristics of the real component of signals
miuX = sum(Xw)./N;
sigmaX = sqrt(sum((Xw-miuX).^2)./N);
miuY = sum(Yw)./N;
sigmaY = sqrt(sum((Yw-miuY).^2)./N);

% Compute cross-correlation coefficient
rou = 1/(N*sigmaX*sigmaY)*sum((Xw-miuX).*(Yw-miuY));

Y = rou;

end