close all
clear
clc

ft = 0.5;
pt = pi;

fc = 0.5;
pc = 0;

N = 256;
Fs = 1000;

idxT = (0:N-1)/Fs;
st = cos(2*pi*ft*idxT+pt);
sc = cos(2*pi*fc*idxT+pc);

stFFT = fft(st,N);
scFFT = fft(sc,N);
Tw1 = abs(stFFT);
Cw1 = abs(scFFT);
Tw = real(stFFT);
Cw = real(scFFT);