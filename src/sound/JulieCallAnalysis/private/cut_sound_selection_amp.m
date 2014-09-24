function [ncuts soundCuts spectroCuts to fo] = cut_sound_selection_amp(sound_in, samprate, amp_env, max_min_ind, duration, pl)
% Cuts the calls into sections based on max_min and centers then on peak
% enveloppe.  Returns calls that are "duration" long. Duration is in
% seconds.

% Length of sound cuts - make sure it is an odd number
cutPts = fix(duration*samprate);
cutPts = cutPts + (1-rem(cutPts, 2));
minPts = fix(0.03*samprate);   % ignore cuts that are below 30 ms.
% maxEnd = fix(0.6*samprate);    % Set at 600 ms for neural data
maxEnd = fix(4*samprate);   % Set here at 4 seconds

% find out the number of valid cuts.
nmax = length(max_min_ind);
ncuts_max = (nmax-1)/2;  
ncuts = 0;
for i = 1:ncuts_max
    ind_beg = -max_min_ind((i-1)*2+1);
    ind_end = -max_min_ind((i-1)*2+3);
    if (ind_end - ind_beg < minPts)
        continue;
    end
    if (ind_end > maxEnd)
        break;
    end
    ncuts = ncuts + 1;
end

% Parameters for the Spectrogram
nstd = 6;
fband = 50;
twindow = 1000*nstd/(fband*2.0*pi);           % Window length in ms - 6 times the standard dev of the gaussian window
winLength = fix(twindow*samprate/1000.0);  % Window length in number of points
winLength = fix(winLength/2)*2;            % Enforce even window length
increment = fix(0.001*samprate);           % Sampling rate of spectrogram in number of points - set at 1 kHz
f_high=12000;  % Upper frequency bound to save relevant section of spectrogram 
DBNOISE = 40;   % Floor threshold for spectrogram

% Generate a blank sound to get dimensions of corresponding spectrogram
sound_temp = zeros(1, cutPts);
[s, to, fo, pg] = GaussianSpectrum(sound_temp, increment, winLength, samprate); 
nt = length(to);
findmax = find(fo >= f_high, 1, 'first');
fo = fo(1:findmax);
nf = length(fo);
clear sound_temp

% Make space of the return values
soundCuts = zeros(ncuts, cutPts);
spectroCuts = zeros(ncuts, nf*nt);

% Cut the sounds and calculate the spectrograms. 
ic = 0;
for i = 1:ncuts_max
    ind_beg = -max_min_ind((i-1)*2+1);
    ind_end = -max_min_ind((i-1)*2+3);
    ind_max = max_min_ind((i-1)*2+2);
    if (ind_end - ind_beg < minPts)
        continue;
    end
    if (ind_end > maxEnd)
        break;
    end
    ic = ic + 1;
    
    % Use the enveloppe profile to find the mean time but only within
    % duration of peak.
    ind_begLimit = ind_max - cutPts;
    ind_endLimit = ind_max + cutPts;
    if ind_beg < ind_begLimit
        ind_beg = ind_begLimit;
    end
    if ind_end > ind_endLimit
        ind_end = ind_endLimit;
    end
    
    sound_temp = amp_env(ind_beg:ind_end);
    ind_mean = 0;
    for j=1:length(sound_temp)
        ind_mean = ind_mean + j*sound_temp(j);
    end
    ind_mean = fix(ind_mean./sum(sound_temp));
    ind_mean = ind_beg - 1 + ind_mean;
    
    ind_beg_tot = ind_mean - (cutPts-1)/2;
    ind_end_tot = ind_mean + (cutPts-1)/2;
    
    if ind_beg_tot < ind_beg
        db = ind_beg - ind_beg_tot + 1;
        ind_beg_final = ind_beg;
    else
        db = 1;
        ind_beg_final = ind_beg_tot;
    end
    if ind_end_tot > ind_end
        de = cutPts - (ind_end_tot-ind_end);
        ind_end_final = ind_end;
    else
        de = cutPts;
        ind_end_final = ind_end_tot;
    end
    
    soundCuts(ic, db:de) = sound_in(ind_beg_final:ind_end_final);
    [s, to, fo, pg] = GaussianSpectrum(soundCuts(ic, :), increment, winLength, samprate); 
    logB = 20*log10(abs(s));
    spectroCuts(ic, :) = reshape(logB(1:nf, :), 1, nf*nt);
    fo = fo(1:nf);
    
    % Plot the spectrogram on figure 2 if the pl flag is set
    if pl        
        figure(2);
        subplot(ncuts,1,ic);
        cla;
        
        maxB = max(max(logB(1:nf,:)));
        minB = maxB-DBNOISE;
        imagesc(to,fo,logB(1:nf,:));          % to is in seconds
        axis xy;
        caxis('manual');
        caxis([minB maxB]);
        
        axis([0 to(end) 0 f_high]);
        ylabel('Frequency kHz');
        cmap = spec_cmap();
        colormap(cmap);
    end
end

return
          

