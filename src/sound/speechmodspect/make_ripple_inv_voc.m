% Makes band-limited noise ripples with uniform distribution or with a distribution given 
% by a modulation sectrum

% Global variables
% User Input/Output variables:
nsongs = 10;         % Number of noise ripples to be made
save_dir = 'D:/Frederic Theunissen/Data/Yale/vocs/syn0_10_1_2';
debug_plot = 1;

% Calculation variables:
max_wx = 2.0;       % Upper cuttoff for spectral frequencies in cycles/KHz
min_wx = 1.0;         % Lower cuttoff for spectral frequencies
max_wt = 10;        % Upper cuttoff for temporal frequencies in Hz
min_wt = 0;        % Lower cuttoff for temporal frequencies
no_it = 30;       % Number of iterations for the spectrographic inversion
f_width = 32;      % Width of the gaussian filters (std) for the spectrographic representation in Hz.  
stim_dur=0.5;         % stimulus duration in s 
samp_rate=50000;    % sampling rate for sound pressure waveform 
amp_samp_rate=1000; % sampling rate for spectrogram
f_low=40;          % Lower bound of filter for sound  
f_high=10000;        % Upper bound of filter for sound
tone_ramp=25;       % tone ramp in ms 
no_of_ripples=100;  % number of componnent ripples in a single noise ripple of stim_dur
log_flg = 1;       % Set to 1 for log spectrum and 0 for linear
mean_amp_voc = 0.549;  % mean amplitude of spectrogram enveloppes in log scale from primate_mod_spect
var_amp_voc = 0.446;   % mean variance in spectrogran enveloppes in log scale from primate_mod_spect
dbatt = 0.0;            % dbattenuation from using entire dynamic range
do_old = 0;             % Set to 1 to also generate sound using the sum of ripple methods

% Calcuate the window length for the spectrogram and the frequency bounds
% of the fft
winLength = fix((1/f_width)*6*samp_rate);    % The window is given by 6 times the std of the gaussian and it is in units of points
if ( rem(winLength, 2) )                % Invert and add assumes an even winLenght
    winLength = winLength +1;
end

fftLen = 0;     
fftLen = max(winLength, fftLen);        % This is the fftLen that will be calculated in ComplexSpectrum
f_step = samp_rate/fftLen;
f_val = 0.0:f_step:samp_rate/2;
nbands = fftLen/2 + 1;                  % This is the number of frequency bands - also equal to the lenght of f_val
increment = samp_rate/amp_samp_rate;
frameCount = floor((stim_dur*samp_rate-winLength)/increment)+1;
t_val = 0:frameCount-1;
t_val = t_val + (winLength*1000.0)/(2*samp_rate);  % This is the time in ms that goes with the spectrogram

% Finds the boundaries in the spectrogram with non-zero entries
for ind_low=1:nbands
    if (f_val(ind_low) >= f_low)
        break;
    end
end
for ind_high=1:nbands
    if (f_val(ind_high) > f_high)
        break;
    end
end
if (ind_high == nbands & f_val(nbands) <= f_high)
    ind_high = nbands;
else
    ind_high = ind_high -1;
end

%memory allocation for calculating the ensemble of song ripples 
mean_amp=zeros(1,nbands);
var_amp=zeros(1,nbands);
maxrip=zeros(1,nsongs);
maxrip_old=zeros(1,nsongs);

%memory allocation for omega_x's and omega_t's selected
omx=zeros(nsongs,no_of_ripples);
omt=zeros(nsongs,no_of_ripples);
omxrip=zeros(nsongs,no_of_ripples);
omtrip=zeros(nsongs,no_of_ripples);


% Generate the unnormalized desired spectrogram
for ns=1:nsongs 
    amp_env_sum = zeros(nbands,frameCount);   
    i=1; 
    
    % Pick random values of spectral and temporal modulations
    while (i <= no_of_ripples) 
        omega_x(i)=rand*(max_wx - min_wx) + min_wx;
        omega_t(i)=rand*(max_wt - min_wt) + min_wt;
        if ( rand < 0.5 )
            omega_t(i) = -omega_t(i);
        end
        % Check for correct bounds: must be below Nyquist sampling rate for temp and spect modulation
        % and within the omega_t*omega_x < i region.  omega_x is in cycles /kHz thus 1000
        if (  abs(omega_x(i)) > 1000./(f_width*2)  | ...
                abs(omega_t(i)) > pi*f_width  ) 
            fprintf(1,'Error: omega_t or omega_x are outside sampled areas\n');
            continue; 
        end
        i=i+1; 
    end 
    
    fprintf(1,'Done finding frequencies for song %d\n', ns);    
    fprintf(1,' Making ripple:\d');
    for i=1:no_of_ripples      
        omx(ns,i)=omega_x(i);
        omt(ns,i)=omega_t(i);
        wrip_x = omega_x(i);
        wrip_t = omega_t(i);
        amp_env = make_amp_env(wrip_x,wrip_t, f_val, frameCount, amp_samp_rate, ind_low, ind_high, 1);
        amp_env_sum = amp_env_sum + amp_env;
        fprintf(1,' %d', i);
        if (mod(i,10) == 0 )
            fprintf(1,'\n   ');
        end
    end    
    mean_amp = mean_amp + mean(amp_env_sum,2)';
    save_file = sprintf('%s/amp_env%d.mat', save_dir, ns);
    save(save_file,'amp_env_sum');
    clear amp_env_sum;
    fprintf(1,'Done making spectrogram for song %d\n\n', ns);
end


% Subtract the mean of envelopes and calculate the variance
mean_amp = mean_amp./nsongs;
mean_mean_amp = mean(mean_amp(ind_low:ind_high)); % This is the overall mean in region of interest

% Make the songs have zero mean and calculate variance
for ns=1:nsongs
   save_file = sprintf('%s/amp_env%d.mat', save_dir, ns);
   load(save_file);
   amp_env_sum=amp_env_sum - mean_mean_amp;
   for i=ind_low:ind_high
		var_amp(i) = var_amp(i) + sum(amp_env_sum(i,:).*amp_env_sum(i,:));
	end
    save(save_file,'amp_env_sum');
    clear amp_env_sum;
end
var_var_amp = sum(var_amp(ind_low:ind_high));
var_amp = var_amp./(nsongs*stim_dur*amp_samp_rate);
var_var_amp = var_var_amp./(nsongs*stim_dur*amp_samp_rate*(ind_high-ind_low+1));
var_fix = sqrt(var_amp_voc/var_var_amp);
amp_fix = mean_amp_voc; 

% Make ripple sound by spetrographic inversion after normalizing desired spectrogram    
fprintf(1, 'Starting Spectrogram Inversions\n');

riperr = zeros(1,nsongs);
for ns=1:nsongs
    save_file = sprintf('%s/amp_env%d.mat', save_dir, ns);
    load(save_file);
    
    % First normalize the spectrogram to have same mean and variance and
    % desired sound ensemble - this sets the contrast of the ripple sounds
    amp_env_sum(ind_low:ind_high,:)=amp_env_sum(ind_low:ind_high,:)*var_fix + amp_fix;      
    
    % This is extra code to prevent negative values for the dc
    if ( log_flg == 0 || log_flg == 1)
        min_env = min(min(amp_env_sum));
        for i=ind_low:ind_high
            amp_env_sum(i,:) = amp_env_sum(i,:) - min_env;
        end
    end
    
    % Rectify amplitude envelope
    [indx,indy]= find(amp_env_sum(ind_low:ind_high,:)< 0.0);
    fprintf(1,'Found %d points below zero for song %d\n', length(indx), ns);
    if(length(indx) > 0)
        i=length(indx);
        amp_env_sum(indx,indy) = 0.0;
    end
    min_env = min(min(amp_env_sum(ind_low:ind_high,:)));
    
    % Transform to exp.
    % Set values outside desired range to min value for better log behavior
    % Not sure this is needed
    set_to_min_flg = 0;
    if ( log_flg )
        if ( set_to_min_flg)
            if ( ind_low ~= 1)
                amp_env_sum(1:ind_low-1, :) =  min_env;
            end
            if ( ind_high ~= nbands)
                amp_env_sum(ind_high+1:nbands, :) =  min_env;
            end
        end
        amp_env_sum = exp(amp_env_sum)-1;
    end
    
    % Save final desired spectrogram as a matlab file
    save(save_file,'amp_env_sum');
    
    % This is a debugging figure
    if ( debug_plot )
        figure(2*ns -1);
        if (do_old)
            subplot(3,1,1);
        else
            subplot(2,1,1);
        end
        if ( log_flg )
            imagesc(t_val, f_val, log(amp_env_sum+1));
        else
            imagesc(t_val, f_val, amp_env_sum);
        end
        axis xy;
        title('Desired Spectrogram');
        
        figure(2*ns);
        if (do_old)
            subplot(3,1,1);
        else
            subplot(2,1,1);
        end
        hold off;
        if ( log_flg )
            plot(t_val,log(amp_env_sum(ind_low+floor((ind_high-ind_low)/2),:)+1));
            hold on;
            plot(t_val, log(amp_env_sum(max(ind_low-5,1),:)+1),'k');
        else
            plot(t_val, amp_env_sum(ind_low+floor((ind_high-ind_low)/2),:));
            hold on;
            plot(t_val, amp_env_sum(max(ind_low-5,1),:),'k');
        end
        xlabel('Time');
        title('Desired Amp enveloppe');
    end
    
    % Random Phase 
    amp_env_phase = (2*rand(size(amp_env_sum))-1) + conj(2*rand(size(amp_env_sum))-1);
    [ripple, theErr] = SpectrumInversion(amp_env_sum, increment, winLength, no_it, amp_env_phase);
    ripple = addramp(ripple, samp_rate, tone_ramp);
    maxr = max(ripple);
    minr = min(ripple);
    maxabs = max(maxr,abs(minr));
    
    ripple_file = sprintf('%s/flatrip%d.mat',save_dir, ns);
    save(ripple_file,'ripple');
    maxrip(ns) = maxabs;
    riperr(ns) = theErr;
    if (do_old)
        fprintf(1,'Starting making ripples by summing sine waves\n');
        ripple_old = make_ripple_sum(amp_env_sum, f_val, t_val, stim_dur, samp_rate, amp_samp_rate);
        ripple_old = addramp(ripple_old, samp_rate, tone_ramp);
        fprintf(1,'Done riplle summation\n');    maxr_old = max(ripple_old);
        minr_old = min(ripple_old);
        maxabs_old = max(maxr_old,abs(minr_old));    
        ripple_file = sprintf('%s/flatrip_old%d.mat',save_dir, ns);
        save(ripple_file,'ripple_old');   
        maxrip_old(ns) = maxabs_old;
    end
    
    for j=1:no_of_ripples
        omxrip(ns,j)=omx(ns,j);
        omtrip(ns,j)=omt(ns,j);
    end
    
    fprintf(1, 'Ripple sound %d finished\n', ns);
    clear amp_env_sum;
end

fprintf(1, 'Mean error for all sounds %g\n', mean(riperr));
% This is the maximum value of all ripples for normalizing
ripmax=max(maxrip);
if (do_old)
    ripmax_old = max(maxrip_old);
end

% Writte ripple sounds as wave file
for ns=1:nsongs
    
    ripple_file= sprintf('%s/flatrip%d.mat',save_dir, ns);
    wav_name = sprintf('%s/flatrip%d.wav',save_dir,ns);
    load(ripple_file);
    ripple = ripple./(ripmax*power(10.0,(dbatt/20))); 
    wavwrite(ripple, samp_rate, 16, wav_name);
    
    if (do_old)
        ripple_file_old = sprintf('%s/flatrip_old%d.mat',save_dir, ns);
        wav_name_old = sprintf('%s/flatrip_old%d.wav',save_dir,ns);
        load(ripple_file_old);
        ripple_old = ripple_old ./(ripmax_old*power(10.0,(dbatt/20))); 
        wavwrite(ripple_old, samp_rate, 16, wav_name_old);
    end
    
    % Debugging figure
    if (debug_plot)    
        figure(2*ns-1);
        amp_ripple = ComplexSpectrum(ripple, samp_rate/amp_samp_rate, winLength);        
        if (do_old) 
            subplot(3,1,2);
        else         
            subplot(2,1,2);
        end
        if (log_flg)
            imagesc(t_val, f_val, log(abs(amp_ripple)+1) );
        else
            imagesc(t_val, f_val, abs(amp_ripple) );
        end
        axis xy;
        title('Spectrogram after inversion');
        
        if (do_old )
            amp_ripple_old = ComplexSpectrum(ripple_old, samp_rate/amp_samp_rate, winLength);
            subplot(3,1,3);
            if (log_flg)
                imagesc(t_val, f_val, log(abs(amp_ripple_old)+1) );
            else
                imagesc(t_val, f_val, abs(amp_ripple_old) );
            end
            axis xy;
            title('Spectrogram after summing ripples - old method');
        end
        
        figure(2*ns);
        if (do_old)
            subplot(3,1,2);
        else
            subplot(2,1,2);
        end
        hold off;
        if (log_flg)
            plot(t_val,log(abs(amp_ripple(ind_low+floor((ind_high-ind_low)/2),:))+1));
            hold on;
            plot(t_val,log(abs(amp_ripple(max(ind_low-5,1),:))+1),'r');
        else
            plot(t_val,abs(amp_ripple(ind_low+floor((ind_high-ind_low)/2),:)));
            hold on;
            plot(t_val,abs(amp_ripple(max(ind_low-5,1),:)),'r');
        end
        xlabel('Time');
        title('Amp enveloppe after inversion');
        
        if (do_old)
            subplot(3,1,3);
            hold off;
            if (log_flg)
                plot(t_val,log(abs(amp_ripple_old(ind_low+floor((ind_high-ind_low)/2),:))+1));
                hold on;
                plot(t_val,log(abs(amp_ripple_old(max(ind_low-5,1),:))+1),'r');
            else
                plot(t_val,abs(amp_ripple_old(ind_low+floor((ind_high-ind_low)/2),:)));
                hold on;
                plot(t_val,abs(amp_ripple_old(max(ind_low-5,1),:)),'r');
            end
            xlabel('Time');
            title('Amp enveloppe after summing');
        end
    end
end