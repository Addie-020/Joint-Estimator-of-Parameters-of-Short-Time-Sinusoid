clear
close all
clc

Fs = 100;                           % Sampling frequency (Hz)
Tt = 2;                             % Total time of sampling (s)
Ns = Tt * Fs;                       % Total sampling points

ft = 0.1;                           % Frequency of test signal (Hz)
pt = pi / 3;                        % Phase of test signal (rad)

xt = (0 : Ns - 1) / Fs;             % Time index
xn = sin(2 * pi * ft * xt + pt);    % Test signal
% Compute mean and variance of test signal
miu0 = sum(xn) / Ns;
sigma0 = sqrt(sum((xn - repmat(miu0, 1, Ns)).^2) / Ns);
% Complement signal information
tn = [xn, miu0, sigma0];

nPop = 200;
X = zeros(2, nPop);
for i = 1 : nPop
    X(1, i) = rand;
    X(2, i) = (rand - 0.5) * pi;
end

tic
Y1 = ObjFun(X, tn, Fs);
time1 = toc;

tic
Y2 = ObjFunComp(X, tn, Fs);
time2 = toc;

time1
time2
