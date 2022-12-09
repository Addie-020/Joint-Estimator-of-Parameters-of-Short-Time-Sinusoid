clear
close all
clc

N = 1024;
Fs = 500;
idxT = (0:N-1)/Fs;
M = 1;
time0 = 0;
iter0 = 0;

for i = 1 : M
    
    fprintf("Iteration No.%d\n", i);
    
    ft = rand(1);
    pt = 2*pi*rand(1);
    at = 1;
    xn = at*cos(2*pi*ft*idxT+pt);
    xn = awgn(xn, 20, 'measured');

    [xBest, K0] = JointEstimator(xn, Fs);
    fe = xBest(1);
    pe = xBest(2);

    errF(i) = fe - ft;
    errP(i) = pe - pt;

    iter0 = iter0 + K0;

end

fmean = mean(abs(errF), 2)';
pmean = mean(abs(errP), 2)';

fprintf('True frequency: %.4f Hz\n', ft);
fprintf('True phase: %.4f (%.4fpi) rad\n', pt, pt/pi);
fprintf('Estimated frequency: %.4f Hz\n', fe);
fprintf('Estimated phase: %.4f (%.4fpi) rad\n', pe, pe/pi);