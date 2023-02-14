close all
clear
clc

Ns = 64;
Fs = 5;
at = 1;
ft = 0.5;
pt = pi/3-0.003;
tIdx = (0:Ns-1)/Fs;
xn = at*cos(2*pi*ft*tIdx+pt);

ac = 1;
fc = 0.5;
pc = pi/3;
% Construct signal
yn = ac*cos(2*pi*fc*tIdx+pc);
% Compute signal means and variances
muX = sum(xn)/Ns;
sigmaX = sqrt(sum((xn-muX).^2)/Ns);
muY = sum(yn)/Ns;
sigmaY = sqrt(sum((yn-muY).^2)/Ns);
% Compute Pearson correlation coefficient
rou = 1/Ns*sum(((xn-muX)/sigmaX).*((yn-muY)/sigmaY));
objVal = 8-exp(rou+1);

rou
objVal