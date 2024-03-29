function [xBest, yBest] = JointEstimatorTraverse(xn, Fs, paramRange, ...
    searchPrec)

%
% Joint estimator of frequency and phase of sinusoid
% Correlation-based method
% Traversing method
% No optimization algorithms
% 
% Input arguments:
%   @xn        : Signal to be estimated
%   @Fs        : Sampling rate
%   @paramRange: Estimation parameter range [fLb, fUb, pLb, pUb]
%   @searchPrec: Searching precision
%
% Output arguments:
%   @xBest  : Optimal point (variable)
%   @fBest  : Optimal value of object function
%
% Author        : Zhiyu Shen @Nanjing University
% Establish Date: Dec 4, 2022
% Revised Data  : Dec 4, 2022
%


%%% Preparation

% Input Vector Size Validation
n = size(xn, 1);
if n ~= 1
    error('Input signal must be in a row vector!');
end

Ns = length(xn);


%%% Search process

% Search range
fLb = paramRange(1);
fUb = paramRange(2);
pLb = paramRange(3);
pUb = paramRange(4);
fInc = searchPrec(1);
pInc = searchPrec(2);

% Iteration
objValMat = CorrSweepTime(fLb, fUb, fInc, pLb, pUb, pInc, xn, Ns, Fs);
[m1, m1Idx] = min(objValMat);
[m2, m2Idx] = min(m1);
iBest = m1Idx(m2Idx)-1;
jBest = m2Idx-1;
yBest = m2;
xBest = [fLb+iBest*fInc, pLb+jBest*pInc];

end % end: function JointEstimator


%%%% Function "CorrSweepTime"

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