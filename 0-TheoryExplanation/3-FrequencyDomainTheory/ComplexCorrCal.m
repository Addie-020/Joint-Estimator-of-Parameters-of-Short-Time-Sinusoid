function Y = ComplexCorrCal(Xw, Yw)
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

% Derive real and image part of signals
reXw = real(Xw);                        % Real component of signal 1
imXw = imag(Xw);                        % Image component of signal 2
reYw = real(Yw);                        % Real component of signal 1
imYw = imag(Yw);                        % Image component of signal 2

% Compute statistical characteristics of the real component of signals
miuX1 = sum(reXw)./N;
sigmaX1 = sqrt(sum((reXw-miuX1).^2)./N);
miuY1 = sum(reYw)./N;
sigmaY1 = sqrt(sum((reYw-miuY1).^2)./N);

% Compute statistical characteristics of the image component of signals
miuX2 = sum(imXw)./N;
sigmaX2 = sqrt(sum((imXw-miuX2).^2)./N);
miuY2 = sum(imYw)./N;
sigmaY2 = sqrt(sum((imYw-miuY2).^2)./N);

% Compute cross-correlation coefficient
rou1 = 1/(N*sigmaX1*sigmaY1)*sum((reXw-miuX1).*(reYw-miuY1));
rou2 = 1/(N*sigmaX2*sigmaY2)*sum((imXw-miuX2).*(imYw-miuY2));

Y = abs(rou1+rou2*1i)/sqrt(2);

end