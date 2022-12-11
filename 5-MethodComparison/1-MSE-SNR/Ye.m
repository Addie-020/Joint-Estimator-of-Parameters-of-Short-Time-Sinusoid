function xBest = Ye(x,Fs,T)

x = x.';
N = length(x);
X = fft(x);
[~,MaxPos] = max(abs(X(2:N/2)).^2);
EstKp = MaxPos;
EstDlt = 0;
EstA = 0;
n = 0:N-1;
for tt = 1:T
    Xp = exp(-1j*2*pi*(EstKp+EstDlt+[-0.5;0.5])/N*n)*x/N;
    Lp = conj(EstA)/N*(1+exp(-1j*4*pi*EstDlt))./(1-exp(-1j*2*pi/N*(2*EstKp+2*EstDlt+[-0.5;0.5])));
    Sp = Xp-Lp;
    EstDlt = EstDlt+N/2/pi*asin(sin(pi/N)*real((Sp(2)+Sp(1))/(Sp(2)-Sp(1))));
    EstA = 1/N*(exp(-1j*2*pi/N*(EstKp+EstDlt)*n)*x-conj(EstA)*(1-exp(-1j*4*pi*EstDlt))/(1-exp(-1j*4*pi/N*(EstKp+EstDlt))));
end
EstP = angle(EstA);
EstF = (EstKp+EstDlt)/N*Fs;

xBest = [EstF, EstP];

end