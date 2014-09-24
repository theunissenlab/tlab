function h = plot_one_spike(spike_time, stim_length, stim_name)

% plots the spectrogram, spike arrival times and average firing rate for
% the spike_time, stim_length and stim_name.  These can be read from
% read_all spikes.  Returns a handle to the figure.

before_stim = 2000;                     % Time before stimulus
total_duration = 6000;                      % Time after stimulus

% The rootcd will be different on the Linux Cluster versus individual
% machines
rootcd = '/auto/fdata/pgill/dataCD_Feb09/all_stims'; % Linux Cluster
%rootcd = 'C:\Users\Frederic\Documents\Data\dataCD_withXF\all_stims'; % Frederic's Paris Laptop

% Get sound, number of trials, etc
ntrials = length(spike_time);
[sound_in, samprate , nbits] = wavread(sprintf('%s/%s', rootcd, stim_name));
if ( ntrials > 50 ) 
	nplot=50;
else
	nplot = ntrials;
end

% Parameters for the Spectrogram
nstd = 6;
fband = 50;
twindow = 1000*nstd/(fband*2.0*pi);           % Window length in ms - 6 times the standard dev of the gaussian window
winLength = fix(twindow*samprate/1000.0);  % Window length in number of points
winLength = fix(winLength/2)*2;            % Enforce even window length
increment = fix(0.001*samprate);           % Sampling rate of spectrogram in number of points - set at 1 kHz
f_low=0;                                 % Lower frequency bounds to get average amplitude in spectrogram
f_high=8000;                               % Upper frequency bound to get average amplitude in spectrogram
DBNOISE = 40;                          % dB in Noise for the log compression - values below will be set to zero.


% Get the figure ready
h=figure;
set(h,'PaperPosition',[0.25 2.5 3.75 8]); % This is for hard copy
set(h,'Position', [50 50 1000 500]);

% Calculate and plot the spectrogram   
subplot(3,1,1);    
[s, to, fo, pg] = GaussianSpectrum(sound_in, increment, winLength, samprate); 
logB = 20*log10(abs(s));
maxB = max(max(logB));
minB = maxB-DBNOISE;            

imagesc(to,fo,logB);          % to is in seconds
axis xy;
caxis('manual');
caxis([minB maxB]); 

after_stim = total_duration - before_stim - stim_length;
v_axis(1) = -before_stim/1000.0;
v_axis(2) = (stim_length + after_stim)/1000.0;
v_axis(3)=f_low; 
v_axis(4)=f_high;

axis(v_axis);                                
xlabel('time (ms)'), ylabel('Frequency');
cmap = spec_cmap();
colormap(cmap);

% plot spike array and calculate average
ntimebins = before_stim + stim_length + after_stim; % Number of time bins in ms
psth = zeros(1, ntimebins);
t=1:ntimebins;
t = (t-before_stim)./1000.0;

for it=1:ntrials
    spike_time_trials = spike_time{it};
    ns = length(spike_time_trials);
    spike_array = zeros(1, ntimebins);
    
    for is=1:ns
        time_ind = round(spike_time_trials(is)) + before_stim;
        if (time_ind < 1 || time_ind > ntimebins)
            fprintf(1, 'Warning time index in plot_one_spike out of bounds: time_ind = %d\n', time_ind);
            continue;
        end
        spike_array(time_ind) = spike_array(time_ind) +1;
    end
    psth = psth + spike_array;
    
    if (it < nplot)
        subplot('position',[0.13 0.66-it*(0.3/nplot) 0.775 (1-0.05)*0.3/nplot]);
        hold on;
        for j=1:ntimebins
            if (spike_array(j) > 0 )
                plot([t(j) t(j)],[0 spike_array(j)],'k');
            end
        end
        axis([-before_stim/1000.0 (stim_length + after_stim)/1000.0 0 1]);
        axis off;
    end
end
psth = psth./ntrials;

subplot(3,1,3)
wind1 = hanning(31)/sum(hanning(31));   % 31 ms smoothing
smpsth = conv(psth,wind1);
plot(t,smpsth(16:length(smpsth)-15)*1000);
axis([-before_stim/1000.0 (stim_length + after_stim)/1000.0 0 1000*max(smpsth)]);
ylabel('Rate (spikes/s)')
xlabel('Time (s)')
