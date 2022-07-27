% Frequency and Phase Co-Estimation

function [f, p, K] = fpestimator(a1, fs)

% 
% Based on correlation
% 
% Adopt conjugate gardient method as 2D research algorithm
% Adopt golden section method as step optimization algorithm
% 

% #codegen

%% Step1: Rough Estimation
% Data pre-processing
N = length(a1);             % Get input sequence length
n = 0 : N - 1;              % Set sequence index
NN = 10 * N;
c = (fft(a1, NN));
cAmp = abs(c);

[maxAmp, uf] = max(cAmp(1:NN/2));
%最高谱峰和次高谱峰

f=f_estimator(a1,fs);
% idx=find(cAmp(1:NN/2)>maxAmp*0.7);
p=0;
if(uf==1)
    f1=(uf-1)*fs/NN;
    f2=uf*fs/NN;
    idx1=uf;
    idx2=uf+1;
else
    if cAmp(uf-1)>cAmp(uf+1)
        f1=(uf-2)*fs/NN;
        f2=(uf-1)*fs/NN;
        idx1=uf-1;
        idx2=uf;
    else
        f1=(uf-1)*fs/NN;
        f2=uf*fs/NN;
        idx1=uf;
        idx2=uf+1;
    end
end
p1=angle(c(idx1));
% p1=correct_p_range(p1);
p2=angle(c(idx2));
% p2=correct_p_range(p2);
if(abs(p1-p2)>5)
    if(p1>p2)
        p2=p2+2*pi;
        p=angle(c(idx2))+(f-f2)*(p1-p2)/(f1-f2);
    else
        p1=p1+2*pi;
        p=angle(c(idx2))+(f-f2)*(p1-p2)/(f1-f2);
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
goal0=target(x0,a1,fs);

xf=[x0(1)+hf,x0(2)];
fd=(target(xf,a1,fs)-goal0)/hf;

xp=[x0(1),x0(2)+hp];
pd=(target(xp,a1,fs)-goal0)/hp;

g0=[fd,pd];
K=0;
f=x0(1);p=x0(2);
while (abs(fd)>epsilon || abs(pd)>epsilon )&& K<20
%     K  
        [x1,goal1,alpha]=updatealpha(x0,g0,a1,goal0,fs);

%  x1
%  goal1

    f=x1(1);p=x1(2);
    xf=[f+hf,p];
    fd=(target(xf,a1,fs)-goal1)/hf;
    xp=[f,p+hp];
    pd=(target(xp,a1,fs)-goal1)/hp;
    g1=[fd,pd];   %求梯度向量g
    r=norm(g1)/norm(g0);
    gc1=g1+r*r*g0;
    if (norm(g1)==0)
        x_est=x1;
        break;
    end

        [x2,goal2,alpha]=updatealpha(x1,gc1,a1,goal1,fs);
%     x2
%     goal2
    f=x2(1);p=x2(2);  
    xf=[f+hf,p];
    fd=(target(xf,a1,fs)-goal2)/hf;

    xp=[f,p+hp];
    pd=(target(xp,a1,fs)-goal2)/hp;

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