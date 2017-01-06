function strf = gaborAud(to, fo, tempCent, freqCent, tempBw, freqBw, tempMod, freqMod, phase)
% Generates a Gabor Spectro-Temporal filter with Cosine Phase.
% Expects t in ms, f in Hz, tempMod is Hz, freqMod in cycles/kHz
% Returns an strf that has dimension fo by to
% Input:
%   to: array of times for which to calcuate strf - in ms
%   fo: array of frequencies for which to calculate strf - in Hz
%   tempCent: center time for temporal Gaussian in ms.
%   freqCent: center frequency for spectral Gaussian in Hz.
%   tempBw: temporal bandwidth in ms.
%   freqBw: spectral bandwidht in Hz.
%   tempMod: temporal Modulation in Hz.
%   freqMod: spectral Modulation in cyc/kHz.

    exparg = zeros(length(fo), length(to));
    cosarg = zeros(length(fo), length(to));
    
    for ifreq = 1:length(fo)
        for it = 1:length(to)
            exparg(ifreq, it) = ((to(it)-tempCent)./tempBw).^2 + ((fo(ifreq)-freqCent)./freqBw).^2;
            cosarg(ifreq, it) = (2*pi/1000.0) * ( (tempMod*to(it)) + (freqMod*fo(ifreq)) );
        end
    end
    
    strf = exp( -exparg ) .* cos(cosarg + phase);
    