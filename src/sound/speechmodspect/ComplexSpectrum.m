function [s, pg] = ComplexSpectrum(input, increment, winLength, doFFTShift, fftLen)
% ComplexSpectrum complex spectrum
% 	s = ComplexSpectrum(input, increment, winLength, doFFTShift)
% 	Compute the complex spectrum of an input signal.
%	Each frame is [winLength]-long and
%	starts [increment] samples after previous frame's start.
%	[doFFTShift] (default value true) controls whether the data is shifted so that
%	it's center is at the start and end of the array (so FFT is more cosine phase)
%	Only zero and the positive frequencies are returned.
%   The hamming window in the original code was changed to a gaussian window where
%   winLength is 6 times the standard deviation of the gaussian. Frederic
%   Theunissen April 2003.

% Windows is made 6 times bigger.
nstd = 6;

%%% ComplexSpectrum %%%
% Malcolm Slaney's code
%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
	fftLen = 0;
end

if nargin < 4
        doFFTShift = 1;
end

if size(input, 1) > 1,
	input = input';
end;

inputLength = length(input);

% Enforce even winLength
if rem(winLength, 2) == 1
    winLength = winLength +1;
end

if inputLength < winLength,
	input(winLength) = 0;
	inputLength = winLength;
end;

frameCount = floor((inputLength-winLength)/increment)+1;

%fftLen = max(2^(nextpow2(winLength)+1), fftLen);
%fftLen = max(2^(nextpow2(winLength)), fftLen);
fftLen = max(winLength, fftLen);

%a = .54;
%b = -.46;
%wr = sqrt(increment/winLength);
%phi = pi/winLength;
%ws = 2*wr/sqrt(4*a*a+2*b*b)*(a + b*cos(2*pi*(0:winLength-1)/winLength + phi));

% The hamming window code is replaced by a Gaussian window code - F.
% Theunissen
wx2 = ((1:winLength)-((winLength+1)/2)).^2;
wvar = (winLength/nstd)^2;
ws = exp(-0.5*(wx2./wvar));


s = zeros(fftLen/2+1, frameCount);
pg = zeros(1, frameCount);
for i=1:frameCount
        start = (i-1)*increment + 1;
        last = start + winLength - 1;
        f = zeros(fftLen, 1);
        f(1:winLength) = ws.*input(start:last);
	pg(i) = std(f(1:winLength));
        if doFFTShift
                f = [f(winLength/2+1:winLength) ; ...
                                     zeros(fftLen-winLength, 1) ; ...
                     f(1:winLength/2)];
        end

        specslice = fft(f);
        s(:,i) = specslice(1:(fftLen/2+1));
end

