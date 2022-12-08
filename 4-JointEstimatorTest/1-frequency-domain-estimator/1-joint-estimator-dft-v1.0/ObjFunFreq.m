function Y = ObjFunFreq(var, Ct, Fs)
%
% Computation of objective function value
% Objective function is based cross correlation coefficient
% X is a nParticles*nvars vector, dimension is shown in row
%
% Input arguments:
%   @var: variables (including frequency and phase component)
%   @Ct : Covariance information of test signal's frequency spectrum
%   @Fs : Sampling rate
%
% Output arguments:
%   @y  : Objective function value of input variable
%
% Author         : Zhiyu Shen @Nanjing University
% Estabilish date: Sept 15, 2022
% Revise data    : Dec 6, 2022
%

% Input vector size validation
[nSig, nVar] = size(var);
if nVar ~= 2
    error('X is not of a valid size!')
end % end: if

% Set parameters
nFFT = Ct(end);                                 % Fetch FFT length
nSam = Ct(end-1);                               % Fetch original signal length
Ct = Ct(1:end-2);                               % Fetch real values of Ct
idxT = (0:nSam-1)/Fs;                           % Time index of samples

% Set freuqency and phase vector
varF = var(:,1);                                % nParticles*1
varP = var(:,2);                                % nParticles*1

% Construct estimating signal
sn = cos(2*pi*varF*idxT+varP);                  % nParticles*NFFT

% Add window and perform DFT on the constructed signal
idxWin = 0:1:nFFT-1;
winSig = 0.54-0.46*cos(2*pi*idxWin/nFFT);       % Window signal (Hamming Window)
snWin = [sn, zeros(nSig,nFFT-nSam)].* ...
    repmat(winSig,nSig,1);                      % Zero padding and add window
snFFT = fft(snWin,nFFT,2);

% Compute the frequency spectrum of constructed signal
Sn1 = abs(snFFT./nFFT);
Sn = Sn1(:,1:nFFT/2);
Sn(:,2:end-1) = 2*Sn(:,2:end-1);

% Compute mean and variance of estimating signal
miuS = sum(Sn,2)/nFFT;                          % nParticles*1
sigmaS = sqrt(sum((Sn-miuS).^2,2)/nFFT);        % nParticles*1

% Compute cross-correlation coefficient (Person correlation coefficient)
Ce = (Sn-miuS)./sigmaS;                         % nParticles*nSequence
Rou = Ce*Ct.'/(nFFT-1);                         % nParticles*1

% Compute objective function value
Y = 8-exp(Rou+1);                               % nParticles*1

end