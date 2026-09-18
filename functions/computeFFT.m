function [Frequency,SolutionFFT]=computeFFT(Time,Solution)
% Compute the one-sided spectrum of a time signal via the Fast Fourier Transform

% Sampling
SamplingPeriod=Time(2)-Time(1);
SamplingFrequency=1/SamplingPeriod;
SignalLength=length(Time);

% Pad the signal to the next power of 2 to improve the performance of the FFT
NumPoints=2^nextpow2(SignalLength);

% Dimension along which the signal evolves in time
Dim=find(size(Solution)==SignalLength,1,'last');

% Two-sided spectrum (extract abs/real/imag of the output as needed)
TwoSided=fft(Solution,NumPoints,Dim)/SignalLength;

% One-sided spectrum
Index=repmat({':'},1,ndims(TwoSided));
Index{Dim}=1:NumPoints/2+1;
SolutionFFT=TwoSided(Index{:});
Index{Dim}=2:NumPoints/2;
SolutionFFT(Index{:})=2*SolutionFFT(Index{:});

% Frequency vector
Frequency=0:SamplingFrequency/NumPoints:SamplingFrequency/2;

end
