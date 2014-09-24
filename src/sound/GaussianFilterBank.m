function [s, fc] = GaussianFilterBank(input, winLength, samprate, flow, fhigh)
%
% Gaussian Filter Bank
% 	s = GuassianFilterBank(input, winLength, samprate)
% 	Decomposes the orginal signal into a set of narrowband filters between flow and fhigh.
%   The width of the narrowband filters is such that the window in time has an std of 1/6 winLength.  
%   The code is written this way to match GaussianSpectrum.
%   The routine returns a matrix s that has the same length as input and the same samprate. 
%   The row of s correspond to the different frequency bands.
%   f0 is the center frequency of each band.  



% Enforce even winLength to have a symmetric window
if rem(winLength, 2) == 1
    winLength = winLength +1;
end

% Make input it into a row vector if it isn't
if size(input, 1) > 1,
	input = input';
end;

% Padd the input with zeros 
pinput = zeros(1,length(input)+winLength);
pinput(winLength/2+1:winLength/2+length(input)) = input;
inputLength = length(pinput);

% Frequency bandwidth
nstd = 6;
fband = samprate*nstd/(winLength*2*pi);        % This line could be deleted if we just gave fband as an argument.

% Find the center frequencies:
fftLen = winLength;         % This is here to match to code for Guassian Spectrum
if rem(fftLen, 2)   % winLength is odd
    select = 1:(fftLen+1)/2;
else
    select = 1:fftLen/2+1;
end
fo = (select-1)'*samprate/fftLen;   
fc = fo(fo>=flow & fo<=fhigh);

% Assign space for s
s = zeros(length(fc), length(input));
    
%Digital filtering
finput = fft(pinput);
if rem(inputLength, 2)   % winLength is odd
    select = 1:(inputLength+1)/2;
else
    select = 1:inputLength/2+1;
end
fi = (select-1)'*samprate/inputLength;


for i=1:length(fc)
    
    % The digital filter for that band
    gf = zeros(size(finput));
    for ig = 1:length(fi)
        wx2 = (fi(ig)-fc(i))^2;
        gf(ig) = exp(-0.5*(wx2./fband^2));
        if (ig > 1)
            gf(end-ig+2) = gf(ig);
        end
    end
        
    % The filtered output
    bfinput = real(ifft(gf.*finput));
    
    % Putting in the array
    s(i,:) = bfinput(winLength/2+1:winLength/2+length(input));
    
end


return

