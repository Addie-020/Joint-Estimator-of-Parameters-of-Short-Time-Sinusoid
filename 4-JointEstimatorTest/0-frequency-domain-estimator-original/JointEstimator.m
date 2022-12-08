function [f, p, K] = JointEstimator(xn, Fs)

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
[nCol, nSam] = size(xn);
if nCol ~= 1
    error('Input signal must be in a row vector!');
end

%%% Step1: Preliminary Estimation (DTFT Peak Search)

% Compute frequency spectrum of the test signal
nFFT = 10*nSam;
xnFFT = fft(xn, nFFT);
Xn = abs(xnFFT);

% Find the peak amplitude of the frequency spectrum
% Correct it with previously proposed method
[~, idxF] = max(Xn);
f = FreqEstimator(xn,Fs);
p = 0;
if idxF == 1
    f1 = (idxF-1)*Fs/nFFT;
    f2 = idxF*Fs/nFFT;
    idx1 = idxF;
    idx2 = idxF+1;
else
    if Xn(idxF-1) > Xn(idxF+1)
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
p1 = angle(xnFFT(idx1));
p2 = angle(xnFFT(idx2));
if(abs(p1-p2)>5)
    if(p1>p2)
        p2=p2+2*pi;
        p=angle(xnFFT(idx2))+(f-f2)*(p1-p2)/(f1-f2);
    else
        p1=p1+2*pi;
        p=angle(xnFFT(idx2))+(f-f2)*(p1-p2)/(f1-f2);
    end
else
    p=p2+(f-f2)*(p1-p2)/(f1-f2);
end
p=correct_p_range(p);
% p=angle(c(idx2))+(f-f2)*(p1-p2)/(f1-f2);

epsilon=10^(-2);

hp=10^(-9);%-8
hf=10^(-9);

x0=[f,p];
goal0=target(x0,xn,Fs);

xf=[x0(1)+hf,x0(2)];
fd=(target(xf,xn,Fs)-goal0)/hf;

xp=[x0(1),x0(2)+hp];
pd=(target(xp,xn,Fs)-goal0)/hp;

g0=[fd,pd];
K=0;
f=x0(1);p=x0(2);
while (abs(fd)>epsilon || abs(pd)>epsilon )&& K<20
    %     K
    [x1,goal1,alpha]=updatealpha(x0,g0,xn,goal0,Fs);

    %  x1
    %  goal1

    f=x1(1);p=x1(2);
    xf=[f+hf,p];
    fd=(target(xf,xn,Fs)-goal1)/hf;
    xp=[f,p+hp];
    pd=(target(xp,xn,Fs)-goal1)/hp;
    g1=[fd,pd];   %求梯度向量g
    r=norm(g1)/norm(g0);
    gc1=g1+r*r*g0;
    if (norm(g1)==0)
        x_est=x1;
        break;
    end

    [x2,goal2,alpha]=updatealpha(x1,gc1,xn,goal1,Fs);
    %     x2
    %     goal2
    f=x2(1);p=x2(2);
    xf=[f+hf,p];
    fd=(target(xf,xn,Fs)-goal2)/hf;

    xp=[f,p+hp];
    pd=(target(xp,xn,Fs)-goal2)/hp;

    g2=[fd,pd];   %求梯度向量g
    x0=x2;goal0=goal2;g0=g2;

    K=K+1;
end
p=correct_p_range(p);
end

function y=target(x,a,fs) %#codegen

% x(1):f
% x(2):p
N=length(a);
% fs=1000;
n=0:N-1;
aa=cos(2*pi*x(1)/fs*n+x(2));
c=fft(a,10*N);
cc=fft(aa,10*N);
y=exp(5*abs(min(min(corrcoef(c,cc)))));
end

%% 校正相位
function p=correct_p_range(p)
p=mod(p,2*pi);
%     if(p>2*pi)
%        p=p-2*pi;
%     else if (p<0)
%            p=p+2*pi;
%         end
%     end
end
%% 0.618法一维搜索
function [x,f,t]=updatealpha(x,g,aa,goal,fs) %#codegen

xlast=x;
g=g/norm(g);
f=goal;
a=0;b=0.05;
xa=x;
fa=goal;
xb=x+b*g;
xb(2)=correct_p_range(xb(2));
fb=target(xb,aa,fs);
k=0;

while(b-a>10^(-9))
    t1=a+0.382*(b-a);
    t2=a+0.618*(b-a);
    x1=x+t1*g;
    x1(2)=correct_p_range(x1(2));
    f1=target(x1,aa,fs);
    x2=x+t2*g;
    x2(2)=correct_p_range(x2(2));
    f2=target(x2,aa,fs);

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