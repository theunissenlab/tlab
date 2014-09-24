function [fund, sal, fund2, soundlen] = fundEstimator(soundIn, fs, logS, t, f)
% Estimates the fundamental of a harmonic sound.
% soundIn is the sound pressure waveformlog spectrogram.
% fs is the sampling rate
% LogS is a log spectrogram calculated at an appropriate scale
% The t and f are the time and frequency that correspond to logS in s and
% Hz.  fund is the fundamental in Hz at the same resolution as the
% spectrogram.
global fundGlobal maxFund minFund

% Some user parameters (should be part of the function at some time)
debugFig = 0;               % Set to zero to eliminate figures.
maxFund = 1500;      % Maximum fundamental frequency
minFund = 300;       % Minimum fundamental frequency
lowFc = 200;         % Low frequency cut-off for band-passing the signal prior to auto-correlation.
highFc = 6000;       % High frequency cut-off
minSaliency = 0.5;   % Threshold in the auto-correlation for minimum saliency - returns NaN for pitch values is saliency is below this number


% Band-pass filtering signal prior to auto-correlation
nfilt = 1024;
if (length(soundIn) < 1024)
    fprintf(1, 'Error in fundEstimator: sound too short for bandpass filtering\n');
end
highpass_filter = fir1(nfilt, 2*lowFc/fs, 'high');
soundIn = filtfilt(highpass_filter, 1, soundIn);
lowpass_filter = fir1(nfilt, 2*highFc/fs);
soundIn = filtfilt(lowpass_filter, 1, soundIn);

if debugFig
    figure(9);
    spec_only(soundIn, 50, fs, 50);
end

% Initializations and useful variables
soundLen = length(soundIn);
nt=length(t);
soundRMS = zeros(1,nt);
fund = zeros(1,nt);
fund2 = zeros(1,nt);
sal = zeros(1,nt);

% Calculate the size of the window for the auto-correlation
alpha = 5; % Number of sd in the Gaussian window
winLen = fix((2*alpha/minFund)*fs);    % Length of Gaussian window based on minFund
if (mod(winLen,2) == 0)  % Make a symmetric window
    winLen = winLen+1;
end
w = gausswin(winLen, alpha);
maxlags = 2*ceil((fs/minFund));
sampS = 1./(t(2)-t(1));


% First calculate the rms in each window

for it = 1:nt
    tval = t(it);   % Center of window in time
    tind = fix(tval*fs);  % Center of window in ind
    tstart = tind - (winLen-1)/2;
    tend = tind + (winLen-1)/2;
    
    if tstart < 1
        winstart = 2 - tstart;
        tstart = 1;
    else
        winstart = 1;
    end
    
    if tend > soundLen
        windend = winLen - (tend-soundLen);
        tend = soundLen;
    else
        windend = winLen;
    end
    
    soundWin = soundIn(tstart:tend).*w(winstart:windend)';
    soundRMS(it) = std(soundWin);
    
end

soundRMSMax = max(soundRMS);

% Calculate the auto-correlation in windowed segments
soundlen = 0;
for it = 1:nt
    
    fund(it) = NaN;
    sal(it) = NaN;
    fund2(it) = NaN;
    if (soundRMS(it) < soundRMSMax*0.1)
        continue;
    end
    soundlen = soundlen + 1;
    tval = t(it);   % Center of window in time
    tind = fix(tval*fs);  % Center of window in ind
    tstart = tind - (winLen-1)/2;
    tend = tind + (winLen-1)/2;
    
    if tstart < 1
        winstart = 2 - tstart;
        tstart = 1;
    else
        winstart = 1;
    end
    
    if tend > soundLen
        windend = winLen - (tend-soundLen);
        tend = soundLen;
    else
        windend = winLen;
    end
    
    soundWin = soundIn(tstart:tend).*w(winstart:windend)';
    
    [autoCorr, lags] = xcorr(soundWin, maxlags, 'unbiased');
    ind0 = find(lags == 0);
    
    % find peaks
    [pksCorr, indPeaksCorr] = findpeaks(autoCorr, 'MINPEAKHEIGHT', max(autoCorr)./10);
    
    % Eliminate center peak and all peaks too close to middle
    
    pksCorr(abs(indPeaksCorr-ind0) < fs/maxFund ) = [];
    indPeaksCorr(abs(indPeaksCorr-ind0) < fs/maxFund ) = [];
    
    % Find max peak
    if isempty(pksCorr)
        pitchSaliency = 0.1; 
    else
        indIndMax = find(pksCorr == max(pksCorr), 1, 'first');
        indMax = indPeaksCorr(indIndMax);   
        fundCorrGuess = fs./abs(lags(indMax));
        pitchSaliency = autoCorr(indMax)./autoCorr(ind0);
    end

    sal(it) = pitchSaliency;
    
    if sal(it) < minSaliency
        continue;
    end
    
    % Calculate the envelope of the auto-correlation after rectification
    envCorr = enveloppe_estimator(autoCorr, fs, maxFund, fs);
    [pksEnvCorr,locsEnvCorr] = findpeaks(envCorr, 'MINPEAKHEIGHT', max(envCorr)./10);
    
    % The max peak should be around zero
    indIndEnvMax = find(pksEnvCorr == max(pksEnvCorr), 1, 'first');
      
    % Take the first peak not in the middle
    if indIndEnvMax+1 > length(locsEnvCorr)
        fundCorrAmpGuess = fundCorrGuess;
        indEnvMax = indMax;
    else
        indEnvMax = locsEnvCorr(indIndEnvMax+1);
        fundCorrAmpGuess = fs./lags(indEnvMax);
    end
    
    % Calculate power spectrum and cepstrum
    Y = fft(soundWin, winLen);
    f = fs/2*linspace(0,1,winLen/2+1);
    flow = find(f >= lowFc, 1, 'first');
    fhigh = find(f >= highFc, 1, 'first');
    
    powSound = 20*log10(abs(Y(1:(winLen+1)/2)));    % This is the power spectrum
    powSoundGood = powSound(1:fhigh);
    maxPow = max(powSoundGood);
    powSoundGood = powSoundGood - maxPow;   % Set zero as the peak amplitude
    powSoundGood(powSoundGood < - 60) = -60;    
    
    % Calculate coarse spectral enveloppe
    warning('off', 'all');
    p = polyfit(f(1:fhigh)', powSoundGood', 3);
    powAmp = polyval(p, f(1:fhigh)'); 
    warning('on', 'all');
    
    % Cepstrum
    CY = dct(powSoundGood-powAmp');            
    
    tCY = 2000.*(0:length(CY)-1)./fs;          % Units of Cepstrum in ms
    fCY = 1000./tCY;                           % Corresponding fundamental frequency in Hz.
    flowCY = find(fCY < lowFc, 1, 'first');
    fhighCY = find(fCY < highFc, 1, 'first');
    
    % Find peak of Cepstrum
    indPk = find(CY(fhighCY:flowCY) == max(CY(fhighCY:flowCY)), 1, 'last');
    indPk = fhighCY + indPk -1; 
    fmass = 0;
    mass = 0;
    indTry = indPk;
    while (CY(indTry) > 0)
        fmass = fmass + fCY(indTry)*CY(indTry);
        mass = mass + CY(indTry);
        indTry = indTry + 1;
        if indTry > length(CY)
            break;
        end
    end
    indTry = indPk - 1;
    if (indTry > 0 )
        while (CY(indTry) > 0)
            fmass = fmass + fCY(indTry)*CY(indTry);
            mass = mass + CY(indTry);
            indTry = indTry - 1;
            if indTry < 1
                break;
            end
        end
    end
    fGuess = fmass/mass;
    if (fGuess == 0  || isnan(fGuess) || isinf(fGuess) )              % Failure of cepstral method
        fGuess = fundCorrGuess;
    end
    fundCepGuess = fGuess;
    
    % Force fundamendal to be bounded
    if (fundCepGuess > maxFund )
        i = 2;
        while(fundCepGuess > maxFund)
            fundCepGuess = fGuess/i;
            i = i + 1;
        end
    elseif (fundCepGuess < minFund)
        i = 2;
        while(fundCepGuess < minFund)
            fundCepGuess = fGuess*i;
            i = i + 1;
        end
    end
    
 
    
    % Fit a Guassian harmonic stack
    fundGlobal = fundCepGuess;
    maxPow = max(powSoundGood'-powAmp);
    
    warning('off', 'all');
    fundFitCep = NonLinearModel.fit(f(1:fhigh)', powSoundGood'-powAmp, @synSpect, [fundCepGuess ones(1,9).*log(maxPow)]); 
%     modelPowCep = synSpect(double(fundFitCep.Coefficients(:,1)), f(1:fhigh));  
    modelPowCep = synSpect(fundFitCep.Coefficients.Estimate, f(1:fhigh));
    errCep = sum((powSoundGood - powAmp' - modelPowCep).^2);
    
    fundFitCep2 = NonLinearModel.fit(f(1:fhigh)', powSoundGood'-powAmp, @synSpect, [fundCepGuess*2 ones(1,9).*maxPow]); 
    modelPowCep2 = synSpect(fundFitCep2.Coefficients.Estimate, f(1:fhigh));
    errCep2 = sum((powSoundGood - powAmp' - modelPowCep2).^2);
    
    warning('on', 'all');
    
    if errCep2 < errCep
        fundFitCep = fundFitCep2;
        modelPowCep =  modelPowCep2;
    end
    fundCepGuess = fundFitCep.Coefficients.Estimate(1);
    if (fundCepGuess > maxFund || fundCepGuess < minFund )
        fundCepGuess = NaN;
    end
    
    % A second cepstrum for the second voice
%     CY2 = dct(powSoundGood-powAmp'- modelPowCep);  
                
    fund(it) = fundCepGuess;        
    
    if (~isnan(fundCepGuess))
        powLeft = powSoundGood- powAmp'- modelPowCep;
        maxPow2 = max(powLeft);
        f2 = 0;
        if ( maxPow2 > maxPow*0.5)    % Possible second peak in central area as indicator of second voice.
            f2 = f(find(powLeft == maxPow2, 1));
            if ( f2 > 1000 && f2 < 4000)
                if (pitchSaliency > minSaliency)
                    fund2(it) = f2;
                end
            end
        end
    end
        
        
    
%     modelPowCorrAmp = synSpect(double(fundFitCorrAmp.Coefficients(:,1)), f(1:fhigh));
%     
%     errCorr = sum((powSoundGood - powAmp' - modelPowCorr).^2);
%     errCorrAmp = sum((powSoundGood - powAmp' - modelPowCorrAmp).^2);
%     errCorrSum = sum((powSoundGood - powAmp' - (modelPowCorr+modelPowCorrAmp) ).^2);
%       
%     f1 = double(fundFitCorr.Coefficients(1,1));
%     f2 = double(fundFitCorrAmp.Coefficients(1,1));
%     
%     if (pitchSaliency > minSaliency)
%         if (errCorr < errCorrAmp)
%             fund(it) = f1;
%             if errCorrSum < errCorr
%                 fund2(it) = f2;
%             end
%         else
%             fund(it) = f2;
%             if errCorrSum < errCorrAmp
%                 fund2(it) = f1;
%             end
%         end
%         
%     end
    
    
    
    if (debugFig )
        figure(10);
        subplot(4,1,1)
        plot(soundWin);
%         f1 = double(fundFitCorr.Coefficients(1,1));
%         f2 = double(fundFitCorrAmp.Coefficients(1,1));
        
        title(sprintf('Saliency = %.2f Pitch AC = %.2f (Hz)  Pitch ACA = %.2f Pitch C %.2f (Hz)\n', pitchSaliency, fundCorrGuess, fundCorrAmpGuess, fundCepGuess));
        
        subplot(4,1,2);
        plot(1000.*lags./fs, autoCorr);
        hold on;
        plot([1000.*lags(indMax)./fs 1000.*lags(indMax)./fs], [0 autoCorr(ind0)], 'k');
        plot(1000.*lags./fs, envCorr, 'r', 'LineWidth', 2);
        plot([1000.*lags(indEnvMax)./fs 1000.*lags(indEnvMax)./fs], [0 autoCorr(ind0)], 'g');
        xlabel('Time (ms)');
        hold off;
        
        
        subplot(4,1,3);
        plot(f(1:fhigh),powSoundGood)
        axis([0 highFc -60 0]);
        
        hold on;
        plot(f(1:fhigh), powAmp, 'b--');
        plot(f(1:fhigh), modelPowCep + powAmp', 'k');
        % plot(f(1:fhigh), modelPowCorrAmp + powAmp', 'g');
        
        for ih=1:6
            plot([fundCorrGuess*ih fundCorrGuess*ih], [-60 0], 'r');
            plot([fundCepGuess*ih fundCepGuess*ih], [-60 0], 'k');
        end
        
        if f2 ~= 0 
            plot([f2 f2], [-60 0], 'g');
        end
            
        xlabel('Frequency (Hz)');
        % title(sprintf('Err1 = %.1f Err2 = %.1f', errCorr, errCorrAmp));
        hold off;
        
        subplot(4,1,4);
        plot(tCY, CY);
        hold on;
%         plot(tCY, CY2, 'k--');
        plot([1000/fundCorrGuess 1000/fundCorrGuess], [0 max(CY)], 'r');
        
        plot([1000/fundCepGuess 1000/fundCepGuess], [0 max(CY)], 'k');
        
        %         plot([(pkClosest-1)/fs (pkClosest-1)/fs], [0 max(CY)], 'g');
        %         if ~isempty(ipk2)
        %             plot([(pk2-1)/fs (pk2-1)/fs], [0 max(CY)], 'b');
        %         end
        %         for ip=1:length(pks)
        %             plot([(locs(ip)-1)/fs (locs(ip)-1)/fs], [0 pks(ip)/4], 'r');
        %         end
        axis([0 1000*length(CY)./(2*fs) 0 max(CY)]);
        xlabel('Time (ms)');
        hold off;
        pause();
    end
end

function synS = synSpect(b, x)
% Generates a faction spectrum made out of gaussian peaks
% fund, sigma, pkmax, dbfloor
global fundGlobal maxFund minFund

npeaks = length(b)-1;
% amp = 25;      % Force 25 dB peaks
sdpk = 60;    % Force 80 hz width

synS = zeros(size(x));


for i=1:npeaks
    a = b(i+1);   % To inforce positive peaks only
    synS = synS + a*exp(-(x-b(1)*i).^2/(2*sdpk^2));
end

if (sum(isinf(synS)) + sum(isnan(synS)))
    for i=1:npeaks
       fprintf(1,'%f ', exp(b(i+1)));   
    end
end


% function synS = synSpect(b, x)
% % Generates a faction spectrum made out of gaussian peaks
% % fund, sigma, pkmax, dbfloor
% 
% npeaks = 6;
% amp = 25;      % Force 25 dB peaks
% sdpk = 100;    % Force 100 hz width
% 
% if (b(1) > 1500)
%     b(1) = 1500;
% end
% 
% synS = ones(size(x)).*b(2);
% 
% for i=1:npeaks
%     synS = synS + amp*exp(-(x-b(1)*i).^2/(2*sdpk^2));
% end

% function synS = synSpect(b, x)
% % Generates a faction spectrum made out of gaussian peaks
% % fund, sigma, pkmax, dbfloor
% 
% npeaks = 8;
% synS = ones(size(x)).*(b(3)-b(4));
% 
% b2 = 100;
% for i=1:npeaks
%     synS = synS + b(4)*exp(-(x-b(1)*i).^2/(2*b2^2));
% end



% % Eliminate fundamental values that are isolated in the histogram
% [nelements,xcenters] = hist(fund);
% indMax = find(nelements == max(nelements))
%
% for i=indMax-1:-1:1
%     if nelements(i) == 0
%         fund(fund<xcenters(i)) = NaN;
%         break;
%     end
% end
%
% for i=indMax+1:10
%     if nelements(i) == 0
%         fund(fund>xcenters(i)) = NaN;
%         break;
%     end
% end





