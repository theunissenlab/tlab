function [xvals, yvals, modspec] = strf_modspec(strf, freqInc)
    
    DBNOISE = 80;

    nb = size(strf, 1);
    nt = size(strf, 2);
    ampsrate = 1000;
    
    % calculate the 2D fft
    fabs = fft2(strf);

    % calculate amplitude and phase
    amp_fabs = abs(fabs);
    phase_fabs = angle(fabs);
    amp_fabs_db = 20*log10(amp_fabs);
    max_amp = max(max(amp_fabs_db));
    min_amp = min(min(amp_fabs_db));

    % Calculate the labels for temporal and spectral frequencies in physical
    % units
    % f_step is the separation between frequency bands
    fstep = freqInc;
    for ib=1:ceil((nb+1)/2)
        dwf(ib)= (ib-1)*(1/(fstep*nb));
        if (ib > 1)
            dwf(nb-ib+2)=-dwf(ib);
        end
    end

    % ampsrate is the sampling rate of the amplitude enveloppe
    for it=1:ceil((nt+1)/2)
        dwt(it) = (it-1)*(ampsrate/nt);
        if (it > 1 )
            dwt(nt-it+2) = -dwt(it);
        end
    end

    
    cmap = spec_cmap();
    imagesc(fftshift(dwt), fftshift(dwf).*1000, fftshift(amp_fabs_db));
    caxis('manual');
    caxis([max_amp-DBNOISE max_amp]);
    title('Amplitude Spectrum');
    axis xy;
    axis([-60 60 0 20]);
    colormap(cmap);
    
    modspec = fftshift(amp_fabs_db);    
    xvals = fftshift(dwt);
    yvals = fftshift(dwf).*1000;
    
    