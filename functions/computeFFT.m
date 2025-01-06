function [Frequency,SolutionFFT]=computeFFT(...
         Time,Solution)
         % Compute Fast Fourier Transform

% Time vector
t=Time;

% Sampling period
ts=t(2)-t(1);

% Sampling frequency
fs=1/ts;

% Signal length
L=length(t);

% New input length from the next power of 2 (improve performance of FFT)
n=2^nextpow2(L);

% Frequency vector
f=0:(fs/n):(fs/2);

% Signal in the time domain
ut=Solution;

% Dimension in which to perform FFT
dim=find(size(ut)==length(t),1,'last');

% Signal in the frequecy domain
uf=fft(ut,n,dim);

% Two-sided spectrum (extract then the quantity of interest, ie abs/real/imag...)
P2=uf/L;

% One-sided spectrum
switch dim
  case 1
    P1=P2(1:n/2+1,:);
    P1(2:end-1,:)=2*P1(2:end-1,:);
  case 2
    P1=P2(:,1:n/2+1);
    P1(:,2:end-1)=2*P1(:,2:end-1);
  case 3
    P1=P2(:,:,1:n/2+1);
    P1(:,:,2:end-1)=2*P1(:,:,2:end-1);
end

% Save FFT results
Frequency=f;
SolutionFFT=P1;

end