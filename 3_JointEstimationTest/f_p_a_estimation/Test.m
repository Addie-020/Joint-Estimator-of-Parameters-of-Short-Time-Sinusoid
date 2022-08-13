clear; clc; close all;
ft = 0.5; pt = 0; at = 1;
Fs =100 ; Tt = 2; Ns = Tt*Fs;
xt = (0 : Ns - 1) / Fs; xn = at * sin(2 * pi * ft * xt + pt);
miu0 = sum(xn) / Ns; sigma0 = sqrt(sum((xn - repmat(miu0, 1, Ns)).^2) / Ns);
Ct = (xn - repmat(miu0, 1, Ns)) / sigma0;
X = [0.5 0.5 0.5; 0 0 0; 1 2 0.5];
Y = ObjFun(X, Ct, Fs);