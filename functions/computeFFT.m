function [Frequency,SolutionFFT]=computeFFT(Time,Solution)
% Compute the one-sided spectrum of a signal sampled at uniform time steps

% Sampling frequency and signal length
SamplingFrequency=1/(Time(2)-Time(1));
Length=length(Time);

% Padded length, taken as the next power of two to improve the performance
LengthPadded=2^nextpow2(Length);

% Frequency vector
Frequency=0:(SamplingFrequency/LengthPadded):(SamplingFrequency/2);

% Dimension along which the signal is sampled
Dim=find(size(Solution)==Length,1,'last');

% Two-sided spectrum (extract afterwards the quantity of interest, e.g. abs/real/imag)
SolutionFFT=fft(Solution,LengthPadded,Dim)/Length;

% One-sided spectrum
Half=repmat({':'},1,ndims(SolutionFFT));
Half{Dim}=1:LengthPadded/2+1;
SolutionFFT=SolutionFFT(Half{:});
Half{Dim}=2:LengthPadded/2;
SolutionFFT(Half{:})=2*SolutionFFT(Half{:});

end
