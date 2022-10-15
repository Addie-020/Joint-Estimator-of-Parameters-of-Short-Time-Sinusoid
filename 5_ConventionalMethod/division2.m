function [ao, fo, po]= division2(a0, f0, p0, fs, N, W)

%%% For causal signal

tIdx = 0 : 1/fs : (N-1)/fs;        %starting from second 0;
xR = a0 * cos(2*pi*f0*tIdx + p0);
xI = a0 * sin(2*pi*f0*tIdx + p0);
x = xR + 0*1i*xI;
x = x .* W';
x = x - mean(x);

%%% FFT
xw = fft(x);
axw = abs(xw);

y0 = max(axw(1:N/2));
I=find(axw==y0);
I=I(1);
f0=(I-1)*fs/N;
axw_max=axw(I);

Rxw=real(xw);
Ixw=imag(xw);
pxw=atan(Ixw./Rxw);
%pxw=phase(xw);
pxw=unwrap(pxw);
%subplot(2,1,2)
%plot((0:N/2-1)*fs/N,pxw(1:N/2))
phi0=pxw(I);
Pi=phi0;
while 0<1;
    if Pi>2*pi;
        Pi=Pi-2*pi;
    elseif Pi<0
        Pi=Pi+2*pi;
    else break;
    end
end
Pi;

%%%%%%%%%%%%%%%%
disp(' ');
%%%Ttft for more accurate
f0_L=f0-fs/N;
f0_H=f0+fs/N;
for k=1:10;
  yL=abs(dtft_f(x,f0_L,N,fs));
  yH=abs(dtft_f(x,f0_H,N,fs));
  if yL>=yH;
    f0_H=f0;
  else
    f0_L=f0;
  end
  f0=(f0_H+f0_L)/2;
end
fo=f0;
y0=dtft_f(x,f0,N,fs)*4/N;
ao=abs(y0);
pxw=phase(y0);
pxw=unwrap(pxw);
%subplot(2,1,2)
%plot((0:N/2-1)*fs/N,pxw(1:N/2))
Pi=pxw;
while 0<1;
    if Pi>2*pi;
        Pi=Pi-2*pi;
    elseif Pi<0
        Pi=Pi+2*pi;
    else break;
    end
end
po=Pi;
%toc
end

function y=dtft_f(x,f0,N,fs)
y=0;
for n=0:N-1;
   y=y+x(n+1)*exp(-i*2*pi*f0/fs*n);
end
end

