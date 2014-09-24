function fund = fundamental_estimator(t, f, logS) 
% Estimates the fundamental from a log spectrogram.
% LogS cand be estimate from GaussiamSpectrum or spec.m
% The t and f are the time and frequency that correspond to logS in s and
% Hz.  fund is the fundamental in Hz.

nt=length(t);
fund = zeros(1,nt);
debug_fig = 0;               % Set to zero to eliminate figures.
max_fund = 4000;      % Maximum fundamental frequency
max_fspect = 5000;    % Maximum frequency to estimate cepstrum and for displaying
max_threshold = 10;   % Threshold to use peaks in dB from maximum at that time slice
min_threshold = 20;    % Threshold in db from maxima to minima
min_pvalue = (0.05)^3;    % if maximum posterior is below this value then use NaN
min_pvalue2 = (0.05)^2;

fband = f(2)-f(1);

% Find the index corresponding to max_fspect
ind_fspect = 0;
for ind_fspect=1:length(f)
    if (f(ind_fspect) > max_fspect )
        break;
    end
end

for it=1:nt
    
    % Correlation analysis: find peaks in power spectrum and figure out if
    % first could be fundamental
    
    max_peak = max(logS(1:ind_fspect, it));
    npeaks = 0;
    clear fpeak apeak;
    % Peaks threshold   
    peak_threshold = max_peak - max_threshold;
 
    new_peak_ok = 1;
    last_min = logS(1, it);

    for ifreq=2:ind_fspect-1
        if (logS(ifreq, it)> logS(ifreq-1, it) && logS(ifreq,it) > logS(ifreq+1, it) )
            if (logS(ifreq, it) < peak_threshold || logS(ifreq,it)-last_min < min_threshold )
                continue;
            end
            if ( new_peak_ok )
                npeaks = npeaks +1;
                fpeak(npeaks) = f(ifreq);
                apeak(npeaks) = logS(ifreq, it);
                %last_min = apeak(npeaks);           % reset the value of the local minima
                new_peak_ok = 0;
            else
                if (logS(ifreq, it) > apeak(npeaks))
                    fpeak(npeaks) = f(ifreq);
                    apeak(npeaks) = logS(ifreq, it);
                    %last_min = apeak(npeaks);       % reset the value of the local minima
                end
            end
        end
        if (logS(ifreq, it)< logS(ifreq-1, it) && logS(ifreq, it)< logS(ifreq+1, it) )
            if (npeaks)
                if (logS(ifreq, it) < apeak(npeaks) - min_threshold )
                    if ( new_peak_ok == 0 )
                       new_peak_ok = 1;
                       last_min = logS(ifreq, it);
                    end
                end
            end
            if (logS(ifreq, it) < last_min )
                last_min = logS(ifreq, it);     %last_min is the last local minima
            end
        end
    end
    
    % check out potential fundamental candidates
    if ( npeaks == 0 )   % None of the peaks make the cuttof - choose the max peak
        npeaks = 1;
        fpeak(1) = f(find(logS(:,it)==max_peak,1));
        apeak(1) = max_peak;
    end
    
    fundvals = diff([0 fpeak]);
    fundstart = min(fundvals);
    fundguess = zeros(1,npeaks);
    for ipk=1:npeaks
        multiplier = round(fpeak(ipk)./fundstart);
        fundguess(ipk) = fpeak(ipk)./multiplier;
    end
    meanfundguess = mean(fundguess);
    stdfundguess = std(fundguess);
    if (stdfundguess < fband)
        stdfundguess = fband;
    end
    
    if (debug_fig)
        figure(10);
        subplot(3,1,1);
        plot(f(1:ind_fspect),logS(1:ind_fspect,it));
        xlabel('Frequency Hz');
        paxis = axis;
        paxis(1)=0;
        paxis(2)=max_fspect;
        axis(paxis);
        title(sprintf('Time slice at %f (ms)', t(it)));
        hold on;
        plot(fpeak, apeak, 'r');
        plot([meanfundguess meanfundguess], [0 max_peak], 'k--');
        hold off;
    end

    % Cepstral analysis
    Hs=spectrum.mtm(1.5);
    Hpsd = psd(Hs,logS(1:ind_fspect,it),'Fs',1./fband);   
    [Pxx Pxxc fdump] = pmtm(logS(1:ind_fspect,it), 1.5,[],1./fband,0.66);
    
    
    if (debug_fig)
        subplot(3,1,2);
        plot(Hpsd);
        xlabel('mSeconds');
    end
    nH = length(Hpsd.data);

    % Set fundamental frequencies to be below max_fund Hz
    for ifreq=2:nH
        if ( 1/Hpsd.frequencies(ifreq) < max_fund )
            iHmin = ifreq;
            break;
        end
    end
    % Obtain likelyhood of fundamendal based on amplitude
    psddata = Hpsd.data(iHmin:end);
    npsddata = length(psddata);
    maxval = max(psddata);         % Max power for spectral frequencies
    %conf = mean(Hpsd.ConfInterval(iHmin:end,2)-Hpsd.ConfInterval(iHmin:end,1));
    conf = mean(Pxxc(iHmin:end,2) - Pxxc(iHmin:end,1));

    plikely = ones(1,npsddata);
    for ips = 1:npsddata
        plikely(ips) = exp(-0.5*(psddata(ips)-maxval).^2./conf.^2);
    end
    
    % And the likelyhood based on the correlation analysis
    plikely_fund = ones(1, npsddata);
    if ( meanfundguess < max_fund )
    for ips=1:npsddata
        ipscorr = ips+iHmin-1;
        sfval = 1./Hpsd.frequencies(ipscorr);
        plikely_fund(ips) = exp(-0.5*(sfval-meanfundguess).^2./stdfundguess.^2);
    end
    end

    % Obtain a prior function
    prior = ones(1, npsddata);  
    if ( it > 2 && ~isnan(fund(it-1)) && ~isnan(fund(it-2)) ) 
        if (it == 3)   % start changing the prior on the third time point
            meansf = (fund(2)+fund(1))/2;   % Take the average of first two points
            stdsf = abs((fund(2)-fund(1)))./sqrt(2);
            if (stdsf < fband)
                stdsf = fband;
            end
        elseif (it > 3)
            meansf = 2*fund(it-1)-fund(it-2);  % Take a derivative
            if ( ~isnan(fund(it-3)) )
                stdsf = std(fund(it-3:it-1))./sqrt(3);
            else
                stdsf = abs((fund(2)-fund(1)))./sqrt(2);
            end
            if (stdsf < fband)
                stdsf = fband;
            end
        end

        for ips=1:npsddata
            ipscorr = ips+iHmin-1;
            sfval = 1./Hpsd.frequencies(ipscorr);
            prior(ips) = exp(-0.5*(sfval-meansf).^2./stdsf.^2);
        end
    end

    post = prior.*plikely.*plikely_fund;
    postmax = max(post);

    if (postmax > min_pvalue )
        imax = find(post==postmax)+iHmin - 1;
        fund(it) = 1.0./Hpsd.frequencies(imax);
    else
        imax = 1;
        fund(it) = NaN;
    end
        
    if (debug_fig)
        hold on;
        plot([Hpsd.frequencies(imax).*1000 Hpsd.frequencies(imax).*1000], [0 10.*log10(maxval)], 'k--');
        hold off;
        psdaxis = axis();

        subplot(3,1,3);
        plot(Hpsd.frequencies(iHmin:end).*1000, plikely, 'b');
        probaxis = axis();
        probaxis(1) = psdaxis(1);
        probaxis(2) = psdaxis(2);
        axis(probaxis);
        hold on;
        plot(Hpsd.frequencies(iHmin:end).*1000, plikely_fund, 'g');
        plot(Hpsd.frequencies(iHmin:end).*1000, prior, 'r');
        plot(Hpsd.frequencies(iHmin:end).*1000, post, 'k');
        legend('Likely Cepstrum', 'Likely Corr', 'Prior', 'Posterior');
        hold off;  
        pause();
    end
 


end

return