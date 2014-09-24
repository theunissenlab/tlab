% Makes band-limited noise ripples with uniform distribution or with a distribution given 
% by a modulation sectrum

% Global variables
% User Input/Output variables:
nsongs = 1;         % Number of noise ripples to be made
ripple_type = 2;    % 1 = song ripple noise 2 = flat ripple noise 
%save_dir = '/home/fet/matlab/stim/flatrip3_100/';
save_dir = 'D:/Frederic Theunissen/Matlab/make_stim/flatrip_limited';
load_dir = 'D:/Frederic Theunissen/Matlab/Amp Mod Power spectra/Data/Birdsong/125';
debug_plot = 1;

% Calculation variables:
max_wx = 2.0;       % Upper cuttoff for spectral frequencies in cycles/KHz
min_wx = 0;         % Lower cuttoff for spectral frequencies
max_wt = 100;        % Upper cuttoff for temporal frequencies in Hz
min_wt = 0;        % Lower cuttoff for temporal frequencies
no_it = 30;       % Number of iterations for the spectrographic inversion
f_width = 400;      % Width of the gaussian filters (std) for the spectrographic representation in Hz.  
stim_dur=2.0;         % stimulus duration in s 
samp_rate=32000;    % sampling rate for sound pressure waveform 
amp_samp_rate=1000; % sampling rate for spectrogram
f_low=250;          % Lower bound of filter for sound  
f_high=8000;        % Upper bound of filter for sound
tone_ramp=25;       % tone ramp in ms 
no_of_ripples=100;  % number of componnent ripples in a single noise ripple of stim_dur
k_corr = 2.57;      % Factor to calculate the effective spectro and temporal bandwidth as explained in Singh and Theunissen
log_flg = 1;       % Set to 1 for log spectrum and 0 for linear

% Calcuate the window length for the spectrogram and the frequency bounds
% of the fft
winLength = (1/f_width)*6*samp_rate;    % The window is given by 6 times the std of the gaussian and it is in units of points
fftLen = 2^(nextpow2(winLength)+1);     % This is the fftLen that will be calculated in ComplexSpectrum
f_step = samp_rate/fftLen;
f_val = 0.0:f_step:samp_rate/2;
nbands = fftLen/2 + 1;                  % This is the number of frequency bands - also equal to the lenght of f_val

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
 
%parameters required for selecting the power 
if (ripple_type == 1)
    power_tol=0.005; 
    cd(load_dir);
    stim_stat_pcolor_with_contour;
    new_stim_stat_fplot=zeros(size(stim_stat_fplot)); 
    %selecting the power from songs power spectrum to calculate ripple 
    max_power = max(max(abs(stim_stat_fplot))); 
    power_cutoff=power_tol*(max_power); 
    
    % this find command give the indices of (x,y) for power values > power_cutoff 
    [w_x,w_t]=find(abs(stim_stat_fplot)> power_cutoff); 
    sum_of_amp=0; 
    for k=1:length(w_x) 
        sum_of_amp=sum_of_amp + (abs(stim_stat_fplot(w_x(k),w_t(k)))); 
    end 
    
    %set up probability density distribution matrix 
    for k=1:length(w_x) 
        new_stim_stat_fplot(w_x(k),w_t(k))= (abs(stim_stat_fplot(w_x(k),w_t(k))))/sum_of_amp; 
    end 
    
    % set up cumulative probability density function 
    cum_power_dist=zeros(size(w_x)); 
    cum_power_dist(1)=new_stim_stat_fplot(w_x(1),w_t(1)); 
    for k=2:length(w_x) 
        cum_power_dist(k)=new_stim_stat_fplot(w_x(k),w_t(k)) + cum_power_dist(k-1); 
    end 
    
end 


% Generate the unnormalized desired spectrogram
for ns=1:nsongs 
    amp_env_sum = zeros(nbands,stim_dur*amp_samp_rate); 
    
    i=1; 
    if (ripple_type == 1)
        while (i <= no_of_ripples) 
            ran_power_val=rand;
            for j=1:size(w_x) 
                if(ran_power_val < cum_power_dist(j)); 
                    omega_x(i)=(w_x(j)-nb+rand-0.5)*(1/(2*fstep*(nb-1)))*1000;% in cycles per KHz 
                    omega_t(i)=(w_t(j)-ntt+rand-0.5)*(1000/(2*(ntt-1))) ;
                    
                    % Check for correct bounds: must be below Nyquist sampling rate for temp and spect modulation
                    % and within the omega_t*omega_x < i region.  omega_x is in cycles /kHz thus 1000
                    
                    if (  abs(omega_x(i)) > 1000./(f_width*k_corr)  | ...
                            abs(omega_t(i)) > k_corr*f_width  ) 
                        fprintf(1,'Error: omega_t or omega_x are outside sampled areas\n');
                        
                        i = i +0;
                    else                    
                        i=i+1;  
                    end
                    break;     
                end 
            end 
        end 
    elseif(ripple_type == 2)
        % Pick random values of spectral and temporal modulations
        while (i <= no_of_ripples) 
            omega_x(i)=rand*(max_wx - min_wx) + min_wx;
            omega_t(i)=rand*(max_wt - min_wt) + min_wt;
            if ( rand < 0.5 )
                omega_t(i) = -omega_t(i);
            end
            % Check for correct bounds: must be below Nyquist sampling rate for temp and spect modulation
            % and within the omega_t*omega_x < i region.  omega_x is in cycles /kHz thus 1000
            if (  abs(omega_x(i)) > 1000./(f_width*k_corr)  | ...
                    abs(omega_t(i)) > k_corr*f_width  ) 
                fprintf(1,'Error: omega_t or omega_x are outside sampled areas\n');
                continue; 
            end
            i=i+1; 
        end 
    end    
    fprintf(1,'Done finding frequencies for song %d\n', ns);
    
    fprintf(1,' Making ripple:\d');
    for i=1:no_of_ripples      
        omx(ns,i)=omega_x(i);
        omt(ns,i)=omega_t(i);
        wrip_x = omega_x(i);
        wrip_t = omega_t(i);
        amp_env = make_amp_env(wrip_x,wrip_t, f_val, stim_dur, amp_samp_rate, ind_low, ind_high, 1);
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
    
% Load the song mean enveloppe and variance

load_file = sprintf('%s/stim_init_avg.avg', load_dir);
load(load_file);
load_file = sprintf('%s/stim_var.avg', load_dir);
load(load_file);
load_file = sprintf('%s/stim_init_count.avg', load_dir);
load(load_file);

flow = stim_init_count(4);
fhigh = stim_init_count(5);
fstep = stim_init_count(6);
f1=flow+fstep/2:fstep:fhigh;
amp_song=interp1(f1,stim_init_avg,f_val(ind_low:ind_high),'nearest', 'extrap');
var_song=interp1(f1,stim_var,f_val(ind_low:ind_high),'nearest', 'extrap');


% Make ripple sound by spetrographic inversion after normalizing desired spectrogram    
fprintf(1, 'Starting Spectrogram Inversions\n');
countrip = 0;
for ns=1:nsongs
    save_file = sprintf('%s/amp_env%d.mat', save_dir, ns);
    load(save_file);
    
    % First normalize the spectrogram to have same mean and variance and
    % desired sound ensemble - this sets the contrast of the ripple sounds
    if ripple_type == 1
        var_fix = sqrt(var_song./var_var_amp);     
        for i=ind_low:ind_high
            amp_env_sum(i,:)=amp_env_sum(i,:)*var_fix(i-ind_low+1) + amp_song(i-ind_low+1);
        end      
    else
        var_fix = sqrt(mean(var_song)/var_var_amp);
        amp_fix = mean(amp_song);
        amp_fix = 0;
        var_fix = sqrt(1.0/var_var_amp);
        amp_env_sum(ind_low:ind_high,:)=amp_env_sum(ind_low:ind_high,:)*var_fix + amp_fix;      
    end    
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
        subplot(3,1,1);
        if ( log_flg )
            imagesc(1:amp_samp_rate*stim_dur, f_val, log(amp_env_sum+1));
        else
            imagesc(1:amp_samp_rate*stim_dur, f_val, amp_env_sum);
        end
        axis xy;
        title('Desired Spectrogram');
        
        figure(2*ns);
        subplot(3,1,1);
        hold off;
        if ( log_flg )
            plot(log(amp_env_sum(ind_low+floor((ind_high-ind_low)/2),:)+1));
            hold on;
            plot(log(amp_env_sum(max(ind_low-5,1),:)+1),'k');
        else
            plot(amp_env_sum(ind_low+floor((ind_high-ind_low)/2),:));
            hold on;
            plot(amp_env_sum(max(ind_low-5,1),:),'k');
        end
        xlabel('Time');
        title('Desired Amp enveloppe');
    end
    
    % Random Phase 
    amp_env_phase = (2*rand(size(amp_env_sum))-1) + conj(2*rand(size(amp_env_sum))-1);
    ripple = SpectrumInversion(amp_env_sum, samp_rate/amp_samp_rate, winLength, no_it, amp_env_phase);
    ripple = addramp(ripple, samp_rate, tone_ramp);
    
    fprintf(1,'Starting making ripples by summing sine waves\n');
    ripple_old = make_ripple_sum(amp_env_sum, f_val, stim_dur, samp_rate, amp_samp_rate);
    ripple_old = addramp(ripple_old, samp_rate, tone_ramp);
    fprintf(1,'Done riplle summation\n');
    %ripple_old = ripple;
    
    maxr = max(ripple);
    minr = min(ripple);
    maxabs = max(maxr,abs(minr));
    
    maxr_old = max(ripple_old);
    minr_old = min(ripple_old);
    maxabs_old = max(maxr_old,abs(minr_old));
    
    if((ripple_type == 1) & (maxabs < 2.5e+10))
        countrip=countrip+1;
        ripple_file =sprintf('%s/songrip%d.mat', save_dir, countrip);
        save(ripple_file,'ripple');
        maxrip(countrip) = maxabs;
        for j=1:no_of_ripples
            omxrip(countrip,j)=omx(ns,j);
            omtrip(countrip,j)=omt(ns,j);
        end
    elseif(ripple_type == 2)
        countrip=countrip+1;
        ripple_file = sprintf('%s/flatrip%d.mat',save_dir, countrip);
        save(ripple_file,'ripple');
        ripple_file = sprintf('%s/flatrip_old%d.mat',save_dir, countrip);
        save(ripple_file,'ripple_old');
        maxrip(countrip) = maxabs;
        maxrip_old(countrip) = maxabs_old;
        for j=1:no_of_ripples
            omxrip(countrip,j)=omx(ns,j);
            omtrip(countrip,j)=omt(ns,j);
        end
    end
    fprintf(1, 'Ripple sound %d finished\n', ns);
    clear amp_env_sum;
end

% This is the maximum value of all ripples for normalizing
ripmax=max(maxrip);
ripmax_old = max(maxrip_old);

% This are the frequency bounds that are used for the spectrogram in STRF estimation 
f_step_desired=62.5;
f_desired = f_low+f_step_desired./2:f_step_desired:f_high;
nskip = round(f_step_desired./f_step);
new_bands = length(f_desired);
mean_amp = zeros(new_bands, 1);


% Writte ripple sounds as wave file, dcp files and spectrograms
for ns=1:countrip
    if(ripple_type == 1)
        ripple_file= sprintf('%s/songrip%d.mat',ns);
        dbatt=0;
        wav_name = sprintf('%s/songrip%d.wav',ns);
        dcp_name = sprintf('%s/songrip%d.dcp',ns);
    elseif(ripple_type == 2)
        ripple_file= sprintf('%s/flatrip%d.mat',save_dir, ns);
        ripple_file_old = sprintf('%s/flatrip_old%d.mat',save_dir, ns);
        dbatt=10;
        wav_name = sprintf('%s/flatrip%d.wav',save_dir,ns);
        dcp_name = sprintf('%s/flatrip%d.dcp',save_dir,ns);
    end
    
    load(ripple_file);
    ripple = ripple./(ripmax*power(10.0,(dbatt/20))); 
    load(ripple_file_old);
    ripple_old = ripple_old ./(ripmax_old*power(10.0,(dbatt/20))); 
    wavwrite(ripple, samp_rate, 16, wav_name);
    write_msong(dcp_name,samp_rate,length(ripple),ripple);
    
    % This is the spectrogram
    amp_ripple = ComplexSpectrum(ripple, samp_rate/amp_samp_rate, winLength);
    amp_ripple_old = ComplexSpectrum(ripple_old, samp_rate/amp_samp_rate, winLength);
    
    % Write out spectrogam in format that can be read by dcp_stim_stat
    fname = sprintf('%s/stim_spect%d.dat', save_dir, ns-1);
    fid = fopen(fname,'w','b');
    ind_new_band = 1;
    amp_log = zeros(new_bands, size(amp_ripple,2));
    for ind=ind_low+nskip./2:nskip:ind_high
        if ( log_flg )
            amp_log(ind_new_band,:) = log(abs(amp_ripple(ind,:))+1);
        else
            amp_log(ind_new_band,:) = abs(amp_ripple(ind,:));
        end
        fwrite(fid, amp_log(ind_new_band,:), 'double');
        ind_new_band = ind_new_band +1;
    end
    mean_amp = mean_amp + mean(amp_log,2);
    
    % Debugging figure
    if (debug_plot)    
        figure(2*ns-1);
        subplot(3,1,2);
        if (log_flg)
            imagesc(log(abs(amp_ripple)+1) );
        else
            imagesc(abs(amp_ripple) );
        end
        axis xy;
        title('Spectrogram after inversion');
        
        subplot(3,1,3);
        if (log_flg)
            imagesc(log(abs(amp_ripple_old)+1) );
        else
            imagesc(abs(amp_ripple_old) );
        end
        axis xy;
        title('Spectrogram after summing ripples - old method');
        
        figure(2*ns);
        subplot(3,1,2);
        hold off;
        if (log_flg)
            plot(log(abs(amp_ripple(ind_low+floor((ind_high-ind_low)/2),:))+1));
            hold on;
            plot(log(abs(amp_ripple(max(ind_low-5,1),:))+1),'r');
        else
            plot(abs(amp_ripple(ind_low+floor((ind_high-ind_low)/2),:)));
            hold on;
            plot(abs(amp_ripple(max(ind_low-5,1),:)),'r');
        end
        xlabel('Time');
        title('Amp enveloppe after inversion');
        
        subplot(3,1,3);
        hold off;
        if (log_flg)
            plot(log(abs(amp_ripple_old(ind_low+floor((ind_high-ind_low)/2),:))+1));
            hold on;
            plot(log(abs(amp_ripple_old(max(ind_low-5,1),:))+1),'r');
        else
            plot(abs(amp_ripple_old(ind_low+floor((ind_high-ind_low)/2),:)));
            hold on;
            plot(abs(amp_ripple_old(max(ind_low-5,1),:)),'r');
        end
        xlabel('Time');
        title('Amp enveloppe after summing');
    end
    
end

% Write a mean amp file that can be read by dcp_stim_stat
mean_amp = mean_amp ./countrip;
fname = sprintf('%s/stim_init_avg.avg', save_dir);
fid = fopen(fname, 'w');
for ind=1:new_bands
    fprintf(fid,'%f\n', mean_amp(ind));
end
fclose(fid);