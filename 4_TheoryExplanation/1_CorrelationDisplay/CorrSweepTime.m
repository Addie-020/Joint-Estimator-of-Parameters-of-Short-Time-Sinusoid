function [R, iMax, jMax] = CorrSweepTime(at, fHead, fEnd, fInc, ...
    pHead, pEnd, pInc, idx, fn, Fs)
%
% Sweep frequency and initial phase
% Estimatie frequency and initial phase among the given scale
% 
% Input arguments:
%   @at   : Amplitude of sinusoid
%   @fHead: Start frequency of sweep
%   @fEnd : End frequency of sweep
%   @fInc : Frequency increment of sweep
%   @pHead: Start phase of sweep
%   @pEnd : End phase of sweep
%   @pInc : Phase increment of sweep
%   @idx  : Sequence index
%   @fn   : Sequence to be measured
%   @Fs   : Sampling frequency
%
% Output arguments:
%   @R   : Correlation coefficient
%   @iMax: Maximum i index
%   @jMax: Maximum j index
%
% Author: Zhiyu Shen @Nanjing University
% Date  : Aug 3, 2022
%

iMax = ceil((fEnd - fHead) / fInc);     % Maximum value of frequency index
jMax = ceil((pEnd - pHead) / pInc);     % Maximum value of phase index
f_num = uint8(iMax + 1);                   % Number of frequency sweeping points
p_num = uint8(jMax + 1);                   % Number of phase sweeping points
R = zeros(f_num, p_num);                    % Generate a matrix to store correlation coefficients

% Cross correlation with for loop
for i = 0 : iMax                               % Outer loop: frequency sweeping
    f_e = fHead + i.*fInc;
    w_e = 2*pi * f_e;
    for j = 0 : jMax-1                         % Inner loop: phase sweeping
        p_e = pHead + j*pInc;
        ge = @(x) at * sin(w_e*x + p_e);
        gn = ge(idx/Fs);
        r = corrcoef(fn, gn);
        R(i + 1, j + 1) = r(1, 2)^2;
    end
end

end