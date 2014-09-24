function [signal,theErr, lastSpectrum] = SpectrumInversion(spectrum, frameIncrement,...
                                        winLength, samprate, iterations, ...
                                        initialPhase, hwr, freq_bounds)
% SpectrumInversion inverse a complex spectrum
%	[signal,lastSpectrum] = SpectrumInversion(spectrum,frameIncrement,...
%                                       winLength, samprate, iterations, ...
%                                       initialPhase, hwr);
%	Invert a spectrogram into the original time-domain signal by recovering the phase.
% 	The [frameIncrement] and [winLength] describe the original parameters of the spectrogram.
%	[iterations] is the number of phase-recovery iterations to perform.
% 	InitialPhase is an array of initial phase guesses (perhaps from a similar reconstruction),
%	and [hwr] is a flag indicating whether the signal is positive only
%	(which can be enforced each time through the loop.)
%	Calling: MatchMagnitudes, InvertAndAdd, GaussianSpectrum
if nargin < 5
        iterations = 10;
end

if nargin < 6
        initialPhase = [];
end

if nargin < 7
        hwr = 0;
end

if nargin < 7
    freq_bounds = -1;
end


iter = 0;                               % First time through
cInitialPhase = 2;                      % Linearly weighted cross-correlation
%cInitialPhase = 0;                     % Zero Phase at start for each frame.
cNoCorrelation = 0;                     % Don't do any cross-correlation
cFFTShift = 0;                          % The FFTShift is not needed since fft followed by ifft gives the right result

if ~isempty(initialPhase)
        spectrum = MatchMagnitudes(spectrum, initialPhase);

        signal=real(InvertAndAdd(spectrum, frameIncrement, winLength, ...
                                        iter, cFFTShift, cNoCorrelation));
else
        signal=real(InvertAndAdd(abs(spectrum), frameIncrement, winLength, ...
                                        iter, cFFTShift, cInitialPhase));
end

if max(max(abs(spectrum))) == 0
        fprintf('Spectrum is zero. Returning without iterations.\n');
        lastSpectrum = spectrum;
        return
end

for iter=1:iterations
    %         lastSpectrum = ComplexSpectrum(signal, frameIncrement, winLength, ...
    %                                         cFFTShift, fftLen);
    
    [lastSpectrum, to, fo, pg] = GaussianSpectrum(signal, frameIncrement, winLength, ...
        samprate);
    if freq_bounds ~= -1       
        maxIndx = find(fo >= freq_bounds(2));
        maxIndx = maxIndx(1);
        minIndx = find(fo < freq_bounds(1));
        minIndx = minIndx(end) + 1;
        lastSpectrum = lastSpectrum(minIndx:maxIndx, :);
    end
    
    theErr = sum(sum((abs(spectrum)-abs(lastSpectrum)).^2))/ ...
        sum(sum(abs(spectrum).^2))*100;
    
    % If error is below 0.1% return
    if (theErr < 0.1)
        break;
    end
    if rem(iter, 5) == 0,
        fprintf('.Residue for iteration %d is %g.%%\n', iter, theErr);
    end;

    lastSpectrum = MatchMagnitudes(spectrum, lastSpectrum);
    signal = real(InvertAndAdd(lastSpectrum, frameIncrement, winLength, ...
        iter, cFFTShift, cInitialPhase));    
    if hwr > 0
        signal = max(0,signal);
    end

end

fprintf('Residue error after %d iterations is %f\n', iter, theErr);

