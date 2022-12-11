function xBest = BaiFine(s,Fs)
s = s.';
N = length(s);
n = (0:N-1).';

%% Step 1
S = fft(s);

%% Step 2
[~,MaxPos] = max(abs(S(2:N/2)).^2);
Kl = MaxPos; % Estimate K_l

% Initial interpolators
%% Step 3 and 4
KlPlusOne = mod(Kl+1,N); % (K_l)+1
KlMinusOne = mod(Kl-1,N); % (K_l)-1
X = S([KlMinusOne+1,Kl+1,KlPlusOne+1]);
r_1 = real(X(2)-X(2+1)); % r_+1
r_M1 = real(X(2)-X(2-1)); % r_-1
i_1 = imag(X(2)-X(2+1)); % i_+1
i_M1 = imag(X(2)-X(2-1)); % i_-1
r = r_1^2+r_M1^2;
i = i_1^2+i_M1^2;
if r > i
    y = r_M1/r_1;
    EstF = abs(Fs/2/pi*acos((y*sin(pi/N*(2*Kl+1))*cos(pi/N*(2*Kl-2))+sin(pi/N*(2*Kl-1))*cos(pi/N*(2*Kl+2)))...
        /(r_M1/r_1*sin(pi/N*(2*Kl+1))+sin(pi/N*(2*Kl-1))))); % Estimate frequency
else
    y = i_M1/i_1;
    a = cos(pi/N*(2*Kl-1))+y*cos(pi/N*(2*Kl+1));
    b = cos(pi/N)+1/2*cos(3*pi/N)+1/2*cos(pi/N*(4*Kl+1))+...
        y*(cos(pi/N)+1/2*cos(3*pi/N)+1/2*cos(pi/N*(4*Kl-1)));
    c = cos(pi/N)*(cos(2*pi/N*(Kl+1))+y*cos(2*pi/N*(Kl-1)));
    EstF = abs(Fs/2/pi*acos((b-sqrt(b^2-4*a*c))/2/a)); % Estimate frequency
end

%% Step 5
R = s.'*exp(-1j*2*pi*n*EstF/Fs);

%% Step 6
y = real(R)/imag(R);
e = sin(2*pi*EstF/Fs)/sin(2*pi*EstF*N/Fs)/cos(2*pi*EstF*(N-1)/Fs);
EstP = atan((N*e+1+y*tan(2*pi*EstF*(N-1)/Fs))/(y*(N*e-1)+tan(2*pi*EstF*(N-1)/Fs))); % Estimate phase

%% Step 7
EstA = real(R*sin(2*pi*EstF/Fs)/(N*exp(1j*EstP)*sin(2*pi*EstF/Fs)+sin(2*pi*EstF*N/Fs)*...
    exp(-1j*(EstP+2*pi*EstF/Fs*(N-1))))); % Estimate amplitude

%% Step 8
EstAP = EstA*exp(1j*EstP); % Get the complex amplitude

% Fine interpolater
%% Step 9
x = 0.1;
Sp = exp(-1j*2*pi*(EstF/Fs+[-x;x]/N)*n.')*s-...
    conj(EstAP)*(1-exp(-1j*2*pi*(2*N*EstF/Fs+[-x;x])))./(1-exp(-1j*2*pi*(2*EstF/Fs+[-x;x]/N)));

%% Step 10
Rk = Sp.*exp(1j*pi*(1-1/N)*[-x;x]);
EstDlt1 = real((Rk(1)-Rk(2))*x/(2*cos(pi*x)*EstAP*N-cos(pi*x/N)*(Rk(1)+Rk(2))));

%% Step 11
EstF = EstF+EstDlt1/N*Fs;

xBest = [EstF, EstP];

end


