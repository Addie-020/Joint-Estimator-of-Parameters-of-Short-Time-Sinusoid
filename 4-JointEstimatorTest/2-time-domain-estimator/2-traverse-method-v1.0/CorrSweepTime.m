function objValMat = CorrSweepTime(fLb, fUb, fInc, ...
    pLb, pUb, pInc, xn, Ns, Fs)
%
% Sweep frequency and initial phase
% Estimatie frequency and initial phase among the given scale
% 
% Input arguments:
%   @fHead    : Start frequency of sweep
%   @fEnd     : End frequency of sweep
%   @fInc     : Frequency increment of sweep
%   @pHead    : Start phase of sweep
%   @pEnd     : End phase of sweep
%   @pInc     : Phase increment of sweep
%   @xn       : Signal to be estimated
%   @Ns       : Number of signal samples
%   @Fs       : Sampling frequency
%
% Output arguments:
%   @objValMat: Matrix storing the values of objective function 
%

iMax = round((fUb-fLb)/fInc);            % Maximum value of frequency index
jMax = round((pUb-pLb)/pInc);            % Maximum value of phase index
fNum = round(iMax);                      % Number of frequency sweeping points
pNum = round(jMax);                      % Number of phase sweeping points
objValMat = zeros(fNum, pNum);           % Generate a matrix to store correlation coefficients
tIdx = (0:Ns-1)/Fs;

% Cross correlation with for loop
ac = 1;
for i = 0 : iMax-1                             % Outer loop: frequency sweep
    fc = fLb + i*fInc;
    for j = 0 : jMax-1                         % Inner loop: phase sweep
        pc = pLb + j*pInc;
        % Construct signal
        yn = ac*cos(2*pi*fc*tIdx+pc);
        % Compute signal means and variances
        muX = sum(xn)/Ns;
        sigmaX = sqrt(sum((xn-muX).^2)/Ns);
        muY = sum(yn)/Ns;
        sigmaY = sqrt(sum((yn-muY).^2)/Ns);
        % Compute Pearson correlation coefficient
        rou = 1/Ns*sum(((xn-muX)/sigmaX).*((yn-muY)/sigmaY));
        objValMat(i+1,j+1) = 8-exp(rou+1);
    end
end

end % end: function CorrSweepTime