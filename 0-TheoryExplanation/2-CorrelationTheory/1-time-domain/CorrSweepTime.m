function [R, iMax, jMax] = CorrSweepTime(fHead, fEnd, fInc, ...
    pHead, pEnd, pInc, L, sigMeas, Fs)
%
% Sweep frequency and initial phase
% Estimatie frequency and initial phase among the given scale
% 
% Input arguments:
%   @fHead  : Start frequency of sweep
%   @fEnd   : End frequency of sweep
%   @fInc   : Frequency increment of sweep
%   @pHead  : Start phase of sweep
%   @pEnd   : End phase of sweep
%   @pInc   : Phase increment of sweep
%   @L      : Signal length
%   @sigMeas: Signal to be measured
%   @Fs     : Sampling frequency
%
% Output arguments:
%   @R   : Correlation coefficient
%   @iMax: Maximum i index
%   @jMax: Maximum j index
%
% Author        : Zhiyu Shen @Nanjing University
% Establish Date: Aug 3, 2022
% Revise Date   : Dec 5, 2022
%

iMax = round((fEnd-fHead)/fInc);            % Maximum value of frequency index
jMax = round((pEnd-pHead)/pInc);            % Maximum value of phase index
fNum = round(iMax);                         % Number of frequency sweeping points
pNum = round(jMax);                         % Number of phase sweeping points
R = zeros(fNum, pNum);                      % Generate a matrix to store correlation coefficients
tIdx = (0:L-1)/Fs;

% Cross correlation with for loop
ac = 1;
for i = 0 : iMax-1                             % Outer loop: frequency sweeping
    fc = fHead + i*fInc;
    for j = 0 : jMax-1                         % Inner loop: phase sweeping
        pc = pHead + j*pInc;
        sigCons = ac*cos(2*pi*fc*tIdx+pc);
        r = corrcoef(sigMeas, sigCons);
        R(i+1,j+1) = r(1,2);
%         R(i+1,j+1) = 8 - exp(r(1,2)+1);
    end
end

end