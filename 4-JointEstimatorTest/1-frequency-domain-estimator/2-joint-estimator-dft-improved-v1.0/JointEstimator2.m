function [xBest, yBest] = JointEstimator2(sigTest, Fs, paramRange)

%
% Joint estimator of frequency and phase of sinusoid
% Frequency domain, DFT-based estimation
% Correlation-based method
% Using conjugate gradient method for searching optimization
%
% Input arguments:
%   @xn         : Signal to be estimated
%   @Fs         : Sampling rate
%   @paramRange : Estimation parameter range [fLb, fUb, pLb, pUb]
%
% Output arguments:
%   @xBest  : Optimal point (variable)
%   @numIter: Number of iterations
%
% Author        : Zhiyu Shen @Nanjing University
% Establish Date: Dec 5, 2022
% Revised Data  : Dec 5, 2022
%

%%% Preparation

% Input Vector Size Validation
[nCol, nSam] = size(sigTest);
if nCol ~= 1
    error('Input signal must be in a row vector!');
end

%%% Step1: Preliminary Estimation (DTFT Peak Search)

% Compute frequency spectrum of the test signal
nFFT = 10*nSam;
testFFT = fft(sigTest, nFFT);
Tw = abs(testFFT);

% Find the peak amplitude of the frequency spectrum
% Correct it with previously proposed method
[~, idxF] = max(Tw);
f0 = FreqEstimator(sigTest,Fs);
if idxF == 1
    f1 = (idxF-1)*Fs/nFFT;
    f2 = idxF*Fs/nFFT;
    idx1 = idxF;
    idx2 = idxF+1;
else
    if Tw(idxF-1) > Tw(idxF+1)
        f1 = (idxF-2)*Fs/nFFT;
        f2 = (idxF-1)*Fs/nFFT;
        idx1 = idxF-1;
        idx2 = idxF;
    else
        f1 = (idxF-1)*Fs/nFFT;
        f2 = idxF*Fs/nFFT;
        idx1 = idxF;
        idx2 = idxF+1;
    end
end

% Compute the angle value at peak and correct it
p1 = angle(testFFT(idx1));
p2 = angle(testFFT(idx2));
if abs(p1-p2) > 5
    if p1 > p2
        p2 = p2 + 2*pi;
        p0 = angle(testFFT(idx2)) + (f0-f2)*(p1-p2)/(f1-f2);
    else
        p1 = p1 + 2*pi;
        p0 = angle(testFFT(idx2)) + (f0-f2)*(p1-p2)/(f1-f2);
    end
else
    p0 = p2 + (f0-f2)*(p1-p2)/(f1-f2);
end
p0 = PhaseCorrect(p0);


%%% Step 2: Precise Estimation (CG Algorithm)

% Set search range
fLb = paramRange(1);
fUb = paramRange(2);
fLbRo = max(fLb, f0-0.05);
fUbRo = min(fUb, f0+0.05);
pLbRo = p0-pi/50;
pUbRo = p0+pi/50;
paramRange = [fLbRo, fUbRo, pLbRo, pUbRo];

% Saerch precise search options
options.maxIter = 5;

% Search for the best estimation
[xBest, yBest, ~] = PreciseEstimator(sigTest, Fs, paramRange, ...
    options, [], []);

end % end: function JointEstimator



%%%% Function "ObjFun"

function y = ObjFun(var, xn, Fs)

%
% Computation of objective function value
% Objective function is based cross correlation coefficient
% X is a nParticles*nvars vector, dimension is shown in row
%
% Input arguments:
%   @var: Variables (including frequency and phase component)
%   @xn : Signal to be estimated
%   @Fs : Sampling rate
%
% Output arguments:
%   @y  : Objective function value of input variable
%

% Construct corresponding signal
N = length(xn);
idxT = (0:N-1)/Fs;
sn = cos(2*pi*var(1)*idxT+var(2));

% Compute DFT of both test and constructed signal
Xw = fft(xn, 10*N);
Sw = fft(sn, 10*N);

% Compute correlation coefficient and objective function value
rou = corrcoef(Xw,Sw);
y = exp(5*abs(rou(1,2)));

end % end: function ObjFun



%%%% Function "PhaseCorrect"

function phasOut = PhaseCorrect(phasIn)

%
% Correct the value of phase to [0,2pi]
%
% Input arguments:
%   @phasIn: Input phase value
%
% Output arguments:
%   @phasOut: Output phase value (corrected)
%

phasOut = mod(phasIn, 2*pi);

end % end: function PhaseCorrect



%%%% Function "GoldenSection"

function [xBest, fBest, alpha] = GoldenSection(x0, g0, sigTest, objVal0, Fs)

%
% Golden section method for 1D search optimization
%
% Input arguments:
%   @x0     : Initial variable value
%   @g0     : Initial gradient value
%   @sigTest: Signal to be estimated
%   @objVal0: Initial objective function value
%   @Fs     : Sampling frequency
%
% Output arguments:
%   @xBest: Optimized variable value
%   @fBest: Optimized objective function value
%   @alpha: Optimized 2D search step length
%

% Set initial values and optimization parameters
g0 = g0/norm(g0);
a = 0;
b = 0.05;

% Initial step
fa = objVal0;
xb = x0 + b*g0;
xb(2) = PhaseCorrect(xb(2));
fb = ObjFun(xb, sigTest, Fs);

% Iteration
while b-a > 10^(-9)
    
    % Compute new step length of both ends
    t1 = a + 0.382*(b-a);
    t2 = a + 0.618*(b-a);
    
    % Update variable value at both ends
    x1 = x0 + t1*g0;
    x1(2) = PhaseCorrect(x1(2));
    f1 = ObjFun(x1, sigTest, Fs);
    x2 = x0 + t2*g0;
    x2(2) = PhaseCorrect(x2(2));
    f2 = ObjFun(x2, sigTest, Fs);
    if f1 > f2
        b = t2;
        fb = f2;
    else
        a = t1;
        fa = f1;
    end

end

% Decide output value
if fa > fb
    fBest = fa;
    alpha = a;
    xBest = x0 + a*g0;
else
    fBest = fb;
    alpha = b;
    xBest = x0 + b*g0;
end

end % end: function GoldenSection