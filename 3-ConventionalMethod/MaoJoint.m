function [xBest, K] = MaoJoint(sigTest, Fs)

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
%   @fBest  : Optimal value of object function
%   @info   : Information of the optimization process
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
p0 = 0;
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

% Define optimization parameters
epsilon = 10^(-2);
hp = 10^(-9);
hf = 10^(-9);

% Compute function and gradient value of initial point
% Function value
x0 = [f0, p0];
objVal0 = ObjFun(x0, sigTest, Fs);
% Gradient value
delXf = [x0(1)+hf, x0(2)];
delFunF = (ObjFun(delXf,sigTest,Fs)-objVal0)/hf;
delXp = [x0(1), x0(2)+hp];
delFunP = (ObjFun(delXp,sigTest,Fs)-objVal0)/hp;
g0 = [delFunF, delFunP];

% Iteration
K = 0;
while ((abs(delFunF) > epsilon) || (abs(delFunP)>epsilon)) && (K < 20)
    % Update step length with 0.618 method
    [x1, goal1, ~] = GoldenSection(x0, g0, sigTest, objVal0, Fs);
    % Update variable values
    f0 = x1(1);
    p0 = x1(2);
    % Update gradient values
    delXf = [f0+hf,p0];
    delFunF = (ObjFun(delXf,sigTest,Fs)-goal1)/hf;
    delXp = [f0,p0+hp];
    delFunP=(ObjFun(delXp,sigTest,Fs)-goal1)/hp;
    g1 = [delFunF, delFunP];
    % Update search direction
    r = norm(g1)/norm(g0);
    gc1 = g1 + r*r*g0;
    if norm(g1) == 0
        x0 = x1;
        break;
    end

    [x2, goal2, ~] = GoldenSection(x1, gc1, sigTest, goal1, Fs);

    f0 = x2(1);
    p0 = x2(2);
    delXf = [f0+hf, p0];
    delFunF = (ObjFun(delXf,sigTest,Fs)-goal2)/hf;

    delXp = [f0,p0+hp];
    delFunP = (ObjFun(delXp,sigTest,Fs)-goal2)/hp;

    g2 = [delFunF,delFunP];
    x0 = x2;
    objVal0 = goal2;
    g0 = g2;

    K = K + 1;

end

xBest = x0;
xBest(2) = PhaseCorrect(p0);


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



%%%% Function "FreqEstimator"

function fPeak = FreqEstimator(sigTest, Fs)

%
% Frequency estimator with DFT peak search
% Dichotomous method
% 
% Input arguments:
%   @sigTest: Signal to be estimated
%   @Fs     : Sampling rate
%
% Output arguments:
%   @fPeak: Frequency of the DFT peak
%

%%% Step1: Rough Estimation

% Calculate the DFT of the test signal
nSam = length(sigTest);
idxT = 0:nSam-1; 
nFFT = 10*nSam;
testFFT = fft(sigTest, nFFT);
Tw = abs(testFFT);

% Find the peak of DFT
[~, idxPeak] = max(Tw);
p1 = angle(testFFT(idxPeak));

% Compute the frequency value of the peak and the secondary peak
f1 = (idxPeak-2)*Fs/nFFT;
f2 = idxPeak*Fs/nFFT;


%%% Step2: Precise Estimation

% Compute the correlation coefficient at 
% the left end of the frequency search range
b1 = cos(2*pi*f1*idxT+p1);
B1 = abs(fft(b1,nFFT));
r1Vex = corrcoef(Tw,B1);
r1 = min(r1Vex(1,2));

% Compute the correlation coefficient at 
% the right end of the frequency search range
b2 = cos(2*pi*f2*idxT+p1);
B2 = abs(fft(b2,nFFT));
r2Vex = corrcoef(Tw,B2);
r2 = min(r2Vex(1,2));

% Find the maximum correlation coefficient with dichotomous method
while f2-f1 > 0.000000001
    fm = (f1+f2)/2;
    bm = cos(2*pi*fm*idxT+p1);
    Bm = abs(fft(bm,nFFT));
    rmVex = corrcoef(Tw,Bm);
    rm = min(rmVex(1,2));
    if r2 > r1 
        f1 = fm;
        r1 = rm;
    else
        f2 = fm;
        r2 = rm;
    end
end
fPeak = (f1+f2)/2;

end % end: function FreqEstimator