close all
clear
clc

% Sampling configuration
Fs = 10;
Ns = 64;

% Generate test signal
ft = 0.5;
pt = pi/3;
at = 1;
tIdx = (0:Ns)/Fs;
sigTest = at*cos(2*pi*ft*tIdx+pt);

% Generate constructed signal
fc = 0.7;
pc = pi/3;
ac = 1;
sigCons = ac*cos(2*pi*fc*tIdx+pc);

% Generate noise
SNRdB = 40;
SNRamp = 10.^(SNRdB./20);
sigmaN = at./SNRamp;



