function [y, fo, meanP, Pupper, Plower, stP]=mtpsd_JN(varargin);
%y should just be the fourier-covariance between columns of x.
%fo is the vector of frequencies.
%meanP, Pupper, and Plower are the jackknifed values of the coherence (which is the magnitue of the coherency);
%function [yo, fo, yJ, yupJ, ylowerJ, stJ]=mtcsg(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers);
% Multitaper Cross-Spectral Density, jacknived estimates and errors
%only meanP(:,1,2), Pupper(:,1,2), Plower(:,1,2) is the correct jack-knife.
%These values are the absolute values of the coherency, to get coherence, these values
%must be squared.
%y is the original cross spectrum without jack-knifing or normalizing to get coherency.
% function A=mtcsd(x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers)
% x : input time series
% nFFT = number of points of FFT to calculate (default 1024)
% Fs = sampling frequency (default 2)
% WinLength = length of moving window (default is nFFT)
% nOverlap = overlap between successive windows (default is WinLength/2)
% NW = time bandwidth parameter (e.g. 3 or 4), default 3
% nTapers = number of data tapers kept, default 2*NW -1
%I've changed this program to output the magnitude of coherency, or sqrt of coherence.
% output yo is yo(f)
%
% If x is a multicolumn matrix, each column will be treated as a time
% series and you'll get a matrix of cross-spectra out yo(f, Ch1, Ch2)
% NB they are cross-spectra not coherences. If you want coherences use
% mtcohere

% Original code by Partha Mitra - modified by Ken Harris
% Also containing elements from specgram.m

% default arguments and that
[x,nFFT,Fs,WinLength,nOverlap,NW,Detrend,nTapers] = mtparam(varargin);
winstep = WinLength - nOverlap;

clear varargin; % since that was taking up most of the memory!

nChannels = size(x, 2);
nSamples = size(x,1);

% check for column vector input
if nSamples == 1 
	x = x';
	nSamples = size(x,1);
	nChannels = 1;
end;

% calculate number of FFTChunks per channel
nFFTChunks = round(((nSamples-WinLength)/winstep));
% turn this into time, using the sample frequency
t = winstep*(0:(nFFTChunks-1))'/Fs;

% allocate memory now to avoid nasty surprises later
y=complex(zeros(nFFT, nChannels, nChannels)); % output array for csd
Py=zeros(nFFT, nChannels, nChannels); % output array for psd's
Periodogram = complex(zeros(nFFT, nTapers, nChannels)); % intermediate FFTs
stP = zeros(nFFT, nChannels, nChannels);
varP = zeros(nFFT, nChannels, nChannels);
Temp1 = complex(zeros(nFFT, nTapers)); %Temps are particular psd or csd values for a frequency and taper
Temp2 = complex(zeros(nFFT, nTapers));
Temp3 = complex(zeros(nFFT, nTapers));

JN = complex(zeros(nFFTChunks,nFFT, nChannels, nChannels));  %jackknifed csd
eJ = complex(zeros(nFFT,1));

% calculate Slepian sequences.  Tapers is a matrix of size [WinLength, nTapers]
[Tapers V]=dpss(WinLength,NW,nTapers, 'calc');

% New super duper vectorized alogirthm
% compute tapered periodogram with FFT 
% This involves lots of wrangling with multidimensional arrays.

TaperingArray = repmat(Tapers, [1 1 nChannels]);
for j=1:nFFTChunks
	Segment = x((j-1)*winstep+[1:WinLength], :);
	if (~isempty(Detrend))
		Segment = detrend(Segment, Detrend);
	end;
	SegmentsArray = permute(repmat(Segment, [1 1 nTapers]), [1 3 2]);
	TaperedSegments = TaperingArray .* SegmentsArray;

	Periodogram(:,:,:) = fft(TaperedSegments,nFFT);

	% Now make cross-products of them to fill cross-spectrum matrix
	for Ch1 = 1:nChannels
		for Ch2 = Ch1:nChannels % don't compute cross-spectra twice
			Temp1 = squeeze(Periodogram(:,:,Ch1));
			Temp2 = squeeze(Periodogram(:,:,Ch2));	
			Temp2 = conj(Temp2);
			Temp3 = Temp1 .* Temp2;

            %eJ and eJ2 are the sum over all the tapers.
            eJ=sum(Temp3, 2)/nTapers;
            JN(j,:,Ch1, Ch2) = eJ/nFFTChunks;
			y(:,Ch1, Ch2)= y(:,Ch1,Ch2) + eJ/nFFTChunks;
		end
	end
end 

% now fill other half of matrix with complex conjugate
for Ch1 = 1:nChannels
	for Ch2 = (Ch1+1):nChannels % don't compute cross-spectra twice
		y(:, Ch2, Ch1) = y(:,Ch1,Ch2);
        Py(:, Ch1, Ch2) = abs(y(:,Ch1,Ch2));
	end
end


for j = 1:nFFTChunks
    JN(j,:, :, :) = abs(y - squeeze(JN(j,:, :,:)));
    for Ch1 = 1:nChannels
        for Ch2 = (Ch1+1):nChannels  
           
            % Obtain the pseudo values
            JN(j,:, Ch1, Ch2) = (nFFTChunks*Py(:, Ch1, Ch2)' - squeeze(JN(j,:, Ch1,Ch2)))/(nFFTChunks-1);
        end
    end  
end

meanP=squeeze(mean(JN,1));

for Ch1=1:nChannels
    for Ch2=Ch1:nChannels
        varP(:,Ch1, Ch2) = var(JN(:,:,Ch1, Ch2),1);
    end
end

%upper and lower bounds will be 2 standard deviations away.
stP=sqrt(varP);
size(stP);
size(meanP);

Pupper = meanP + 2*stP;
Plower = meanP - 2*stP;


		
% set up f array
if ~any(any(imag(x)))    % x purely real
	if rem(nFFT,2),    % nfft odd
		select = [1:(nFFT+1)/2];
	else
		select = [1:nFFT/2+1];
	end
    meanP = meanP(select,:,:);
    Pupper = Pupper(select,:,:);
    Plower = Plower(select,:,:);
	y = y(select,:,:);
else
	select = 1:nFFT;
end
	
fo = (select - 1)'*Fs/nFFT;

% we've now done the computation.  the rest of this code is stolen from
% specgram and just deals with the output stage

if nargout == 0
    % take abs, and plot results
    newplot;
    for Ch1=1:nChannels, 
        for Ch2 = 1:nChannels
            subplot(nChannels, nChannels, Ch1 + (Ch2-1)*nChannels);
            plot(f,20*log10(abs(y(:,Ch1,Ch2))+eps));
            grid on;
            if(Ch1==Ch2) 
                ylabel('psd (dB)'); 
            else 
                ylabel('csd (dB)'); 
            end;
            xlabel('Frequency');
        end
    end
end
% Original mtcsd Written by Kenneth D. Harris
% Jackknife estimates over FFT chunks by Anne Hsu and Frederic Theunissen
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
% any comments, or if you make any extensions
% let me know at harris@axon.rutgers.edu