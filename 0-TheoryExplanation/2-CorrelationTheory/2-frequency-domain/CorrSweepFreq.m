function [R, iMax, jMax] = CorrSweepFreq(fHead, fEnd, fInc, ...
    pHead, pEnd, pInc, Xn, nSam, nFFT, Fs)
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
%   @Xn     : Frequency spectrum of the test signal
%   @nSam   : Number of samples of the original signal
%   @nFFT   : Number of FFT points
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
tIdx = (0:nSam-1)/Fs;

% Cross correlation with for loop
ac = 1;
for i = 0 : iMax-1                             % Outer loop: frequency sweeping
    fc = fHead + i*fInc;
    for j = 0 : jMax-1                         % Inner loop: phase sweeping
        pc = pHead + j*pInc;
        sn = ac*cos(2*pi*fc*tIdx+pc);
        % Add window and perform DFT on the constructed signal
        idxWin = 0:1:nFFT-1;
        winSig = 0.54-0.46*cos(2*pi*idxWin/nFFT);
        snWin = [sn, zeros(1,nFFT-nSam)].*winSig;
        snFFT = fft(snWin,nFFT);
        % Compute the frequency spectrum of constructed signal
        Sn1 = abs(snFFT./nFFT);
        Sn = Sn1(1:nFFT/2);
        Sn(2:end-1) = 2*Sn(2:end-1);
        r = corrcoef(Xn, Sn);
        R(i+1,j+1) = r(1,2);
%         R(i+1,j+1) = 8 - exp(r(1,2)+1);
    end
end

end