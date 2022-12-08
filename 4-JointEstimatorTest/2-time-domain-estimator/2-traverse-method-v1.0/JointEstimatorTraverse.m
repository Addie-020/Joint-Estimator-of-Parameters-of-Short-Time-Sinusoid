function [xBest, yBest] = JointEstimatorTraverse(xn, Fs, paramRange, searchPrec)

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


%%% Compute Sequence Information

% Compute mean and variance of test signal
Ns = length(xn);
miu0 = sum(xn) / Ns;
sigma0 = sqrt(sum((xn-miu0).^2) / Ns);

% Compute signal information for correlation computation
Ct = (xn-miu0) ./ sigma0;


%%% Search process

% Search range
fLb = paramRange(1);
fUb = paramRange(2);
pLb = paramRange(3);
pUb = paramRange(4);
fInc = searchPrec(1);
pInc = searchPrec(2);

% Iteration
xIter = zeros(3, 2);
yIter = 3*zeros(3, 1);
for iter = 1 : 3
    R = CorrSweepTime(fLb, fUb, fInc, pLb, pUb, pInc, Ct, Ns, Fs);
    [m1, m1Idx] = min(R);
    [m2, m2Idx] = min(m1);
    iIter = m1Idx(m2Idx)-1;
    jIter = m2Idx-1;
    yIter(iter) = m2;
    xIter(iter,1) = fLb + iIter*fInc;
    xIter(iter,2) = pLb + jIter*pInc;
end % end: for

xBest = sum(xIter)./3;
yBest = sum(yIter)./3;

end % end: function JointEstimator


%%%% Function "CorrSweepTime"

function R = CorrSweepTime(fLb, fUb, fInc, ...
    pLb, pUb, pInc, Ct, Ns, Fs)
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
%   @Ct     : Covariance information of sequence to be estimated
%   @Ns     : Number of signal samples
%   @Fs     : Sampling frequency
%
% Output arguments:
%   @R   : Correlation coefficient
%
% Author        : Zhiyu Shen @Nanjing University
% Establish Date: Aug 3, 2022
% Revise Date   : Dec 5, 2022
%

iMax = round((fUb-fLb)/fInc);            % Maximum value of frequency index
jMax = round((pUb-pLb)/pInc);            % Maximum value of phase index
fNum = round(iMax);                         % Number of frequency sweeping points
pNum = round(jMax);                         % Number of phase sweeping points
R = zeros(fNum, pNum);                      % Generate a matrix to store correlation coefficients
tIdx = (0:Ns-1)/Fs;

% Cross correlation with for loop
ac = 1;
for i = 0 : iMax-1                             % Outer loop: frequency sweeping
    fc = fLb + i*fInc;
    for j = 0 : jMax-1                         % Inner loop: phase sweeping
        pc = pLb + j*pInc;
        sigCons = ac*cos(2*pi*fc*tIdx+pc);
        % Compute mean and variance of estimating signal
        miuS = sum(sigCons,2)/Ns;
        sigmaS = sqrt(sum((sigCons-miuS).^2, 2)/Ns);
        % Compute cross-correlation coefficient (Person correlation coefficient)
        Ce = (sigCons-miuS)./sigmaS;
        Rou = Ce*Ct.'/Ns;
        R(i+1,j+1) = 8-exp(Rou+1);
    end
end

end % end: function CorrSweepTime