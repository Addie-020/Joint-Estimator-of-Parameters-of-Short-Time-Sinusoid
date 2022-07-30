function [fBest, pBest, K] = fpestimator(X, fs)

%
% Based on correlation
%
% Adopt conjugate gardient method as 2D research algorithm
% Adopt golden section method as step optimization algorithm
%


%%% Step1: Rough Estimation

% Data pre-processing
N = length(X);              % Input sequence length
% n = 0 : N - 1;              % Input sequence index
Nfft = 10 * N;              % FFT points
c = fft(X, Nfft);           % FFT for original signal
cAmp = abs(c);              % Signal's magnitude spectrum

% Find the highest point of magnitude spectrum and its freqeuncy index
[~, uf] = max(cAmp(1 : Nfft / 2));

% Estimation based on the calculation of signal spectrum
fHat = f_estimator(X, fs);

% Determine estimation bound according to signal spectrum
% One of the end of bound is the frequency with maximum magnitude
if uf == 1
    f1 = (uf - 1) * fs / Nfft;
    f2 = uf * fs / Nfft;
    idx1 = uf;
    idx2 = uf + 1;
else
    if cAmp(uf - 1) > cAmp(uf + 1)
        f1 = (uf - 2) * fs / Nfft;
        f2 = (uf - 1) * fs / Nfft;
        idx1 = uf - 1;
        idx2 = uf;
    else
        f1 = (uf - 1) * fs / Nfft;
        f2 = uf * fs / Nfft;
        idx1 = uf;
        idx2 = uf + 1;
    end
end

% Calculate phase of both ends of the bound
p1 = angle(c(idx1));
p2 = angle(c(idx2));

% Estimatie initial phase roughly
if abs(p1 - p2) > 5
    if p1 > p2
        phaNew = p2 + (fHat - f2) * (p1 - (p2 + 2 * pi)) / (f1 - f2);
    else
        phaNew = p2 + (fHat - f2) * ((p1 + 2 * pi) - p2) / (f1 - f2);
    end
else
    phaNew = p2 + (fHat - f2) * (p1 - p2) / (f1 - f2);
end
pHat = correct_p_range(phaNew);


%%% Step 2: Main Estimation Process with Optimization Algorithm

epsilon = 10^(-2);

df = 10^(-9);                               % Infinitesimal of frequency variable
dp = 10^(-9);                               % Infinitesimal of phase variable

x0 = [fHat, pHat];                          % Initial point

% Process the initial point
x = x0;
freqValNewTemp = x(1);
phaValNewTemp = x(2);
funVal = objFun(x, X, fs);
% Calculate partial differential on frequency
xIncFreq = [freqValNewTemp + df, phaValNewTemp];
funValIncFreq = objFun(xIncFreq, X, fs);
gFreq = (funValIncFreq - funVal) / df;
% Calculate partial differential on phase
xIncPha = [freqValNewTemp, phaValNewTemp + dp];
funValIncPha = objFun(xIncPha, X, fs);
gPha = (funValIncPha - funVal) / df;
% Complete differential
g = [gFreq, gPha];

K = 0;
while (abs(gFreq)>epsilon || abs(gPha)>epsilon )&& K<20

    [xNewTemp, funValNewTemp, ~] = GoldenSection(x, g, X, funVal, fs);
    freqValNewTemp = xNewTemp(1);
    phaValNewTemp = xNewTemp(2);
    
    % Calculate partial differential on frequency
    xIncFreq = [freqValNewTemp + df, phaValNewTemp];
    funValIncFreq = objFun(xIncFreq, X, fs);
    gFreqNewTemp = (funValIncFreq - funVal) / df;
    % Calculate partial differential on phase
    xIncPha = [freqValNewTemp, phaValNewTemp + dp];
    funValIncPha = objFun(xIncPha, X, fs);
    gPhaNewTemp = (funValIncPha - funVal) / df;
    % Complete differential    
    gNewTemp = [gFreqNewTemp, gPhaNewTemp];

    r = norm(gNewTemp) / norm(g);
    gc1 = gNewTemp + r * r * g;
    
    xBest = xNewTemp;
    if (norm(gNewTemp) == 0)
        break;
    end

    [xNew, funValNew, ~] = GoldenSection(x, gc1, X, funValNewTemp, fs);
    freqNew = xNew(1);
    phaNew = xNew(2);

    % Calculate partial differential on frequency
    xIncFreq = [freqNew + df, phaNew];
    funValIncFreq = objFun(xIncFreq, X, fs);
    gFreqNew = (funValIncFreq - funVal) / df;
    % Calculate partial differential on phase
    xIncPha = [freqNew, phaNew + dp];
    funValIncPha = objFun(xIncPha, X, fs);
    gPhaNew = (funValIncPha - funVal) / df;
    % Complete differential    
    gNew = [gFreqNew, gPhaNew];

    x = xNew;
    funVal = funValNew;
    g = gNew;

    K=K+1;

end

fBest = xBest(1);
pBest = correct_p_range(xBest(2));

end



%% Object Function
function y = objFun(x, a, fs) %#codegen

% x(1):f
% x(2):p

N = length(a);
n = 0 : N - 1;
aa = cos(2 * pi * x(1) / fs * n + x(2));
c = fft(a, 10 * N);
cc = fft(aa, 10 * N);
y = exp(5 * abs(min(min(corrcoef(c, cc)))));

end



%% Phase Correction

function pOut = correct_p_range(pIn)

pOut = mod(pIn, 2 * pi);
%     if(p>2*pi)
%        p=p-2*pi;
%     else if (p<0)
%            p=p+2*pi;
%         end
%     end

end



%% Golden Section Method

function [x, f, t] = GoldenSection(x, g, aa, goal, fs)

xlast = x;
g = g / norm(g);
f = goal;
a = 0;
b = 0.05;
xa = x;
fa = goal;
xb = x + b*g;
xb(2) = correct_p_range(xb(2));
fb = objFun(xb,aa,fs);
k = 0;

while b-a > 10^(-9)
    t1 = a + 0.382 * (b - a);
    t2=a+0.618*(b-a);
    x1=x+t1*g;
    x1(2)=correct_p_range(x1(2));
    f1=objFun(x1,aa,fs);
    x2=x+t2*g;
    x2(2)=correct_p_range(x2(2));
    f2=objFun(x2,aa,fs);

    if(f1>f2)
        b=t2;fb=f2;
    else
        a=t1;fa=f1;
    end
    k=k+1;
end
if(fa>fb)
    f=fa;
    t=a;
    x=xlast+a*g;
else
    f=fb;
    t=b;
    x=xlast+b*g;
end

end