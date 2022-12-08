close all
clear
clc

% Sampling configuration
Fs = 10;
Ns = 8;

% Set frequency and phase range
fLb = 0;
fUb = 1;
pLb = 0;
pUb = 2*pi;

% Generate test signal
ft = 0.5;
pt = pi/3;
at = 1;
tIdx = (0:Ns-1)/Fs;
sigTest = at*cos(2*pi*ft*tIdx+pt);

% Generate constructed signal
fc = fLb + 0.01*randi([0 round(100*(fUb-fLb))]);
pc = pLb + 0.01*randi([0 round(100*(pUb-pLb))]);
ac = 1;
sigCons = ac*cos(2*pi*fc*tIdx+pc);

% Generate noise
SNRdB = 10;
SNRamp = 10.^(SNRdB./20);
sigmaN = at./(sqrt(2)*SNRamp);
sigNois = sigmaN*randn(1,Ns);

sigMeas = sigTest + sigNois;
sigEsti = sigTest;

% Compute correlation coefficient
r1Vec = corrcoef(sigTest, sigNois);
r1 = r1Vec(1,2);
r2Vec = corrcoef(sigTest, sigCons);
r2 = r2Vec(1,2);
r3Vec = corrcoef(sigNois, sigCons);
r3 = r3Vec(1,2);
r4Vec = corrcoef(sigMeas, sigCons);
r4 = r4Vec(1,2);
r5Vec = corrcoef(sigEsti, sigNois);
r5 = r5Vec(1,2);
r6Vec = corrcoef(sigEsti, sigTest);
r6 = r6Vec(1,2);
r7Vec = corrcoef(sigEsti, sigMeas);
r7 = r7Vec(1,2);
fprintf("R(s_t,s_N) = %.4f\n", r1);
fprintf("R(s_t,s_c) = %.4f\n", r2);
fprintf("R(s_N,s_c) = %.4f\n", r3);
fprintf("R(s_m,s_c) = %.4f\n", r4);
fprintf('Good Estimation:\n')
fprintf("R(s_e,s_N) = %.4f\n", r5);
fprintf("R(s_e,s_t) = %.4f\n", r6);
fprintf("R(s_e,s_m) = %.4f\n", r7);

% Compute mean and variance of test signal
Ns = length(sigMeas);
miu0 = sum(sigMeas) / Ns;
sigma0 = sqrt(sum((sigMeas-miu0).^2) / Ns);
% Compute signal information for correlation computation
Ct = (sigMeas-miu0) ./ sigma0;
% Compute mean and variance of estimating signal
miuS = sum(sigEsti,2)/Ns;
sigmaS = sqrt(sum((sigEsti-miuS).^2, 2)/Ns);
% Compute cross-correlation coefficient (Person correlation coefficient)
Ce = (sigEsti-miuS)./sigmaS;
r0 = Ce*Ct.'/Ns;

fprintf("r0 = %.4f\n", r0);

