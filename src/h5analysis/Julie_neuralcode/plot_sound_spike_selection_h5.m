function plot_sound_spike_selection_h5(response,sound_in, samprate, on)
before_stim = 1000;                     % Time before stimulus
total_duration = 4500;                      % Time after stimulus
% Get sound, number of trials, etc
ntrials = length(response.trials);
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
figure(1);
h = gcf();
set(h,'PaperPosition',[0.25 2.5 3.75 8]); % This is for hard copy
set(h,'Position', [50 50 1000 500]);

% Calculate and plot the spectrogram   
subplot(4,1,1); 
cla;
[s, to, fo, pg] = GaussianSpectrum(sound_in, increment, winLength, samprate); 
logB = 20*log10(abs(s));
maxB = max(max(logB));
minB = maxB-DBNOISE;            

imagesc(to,fo,logB);          % to is in seconds
axis xy;
caxis('manual');
caxis([minB maxB]); 

after_stim = total_duration - before_stim - 1000*response.stim_duration;
v_axis(1) = -before_stim/1000.0;
v_axis(2) = (1000*response.stim_duration + after_stim)/1000.0;
v_axis(3)=f_low; 
v_axis(4)=f_high;

axis(v_axis);                                
ylabel('Frequency kHz');
cmap = spec_cmap();
colormap(cmap);



% plot spike array and calculate average
ntimebins = before_stim + 1000*response.stim_duration + after_stim; % Number of time bins in ms
psth = zeros(1, ntimebins);
t=1:ntimebins;
t = (t-before_stim)./1000.0;

for it=1:ntrials
    trial = response.trials{it};
    spike_time_trials = trial.spikeTimes;
    ns = length(spike_time_trials);
    spike_array = zeros(1, ntimebins);
    
    for is=1:ns
        time_ind = round(spike_time_trials(is)*1000) + before_stim;
        if (time_ind < 1 || time_ind > ntimebins)
            fprintf(1, 'Warning time index out of bounds for stim# %s trial %d: time_ind = %d\n', response.number, it, time_ind);
            continue;
        end
        spike_array(time_ind) = spike_array(time_ind) +1;
    end
    psth = psth + spike_array;
    
    if (it <= nplot)
        subplot('position',[0.13 0.66-(it-3)*(0.24/nplot) 0.775 (1-0.05)*0.24/nplot]);
        cla;
        hold on;
        for j=1:ntimebins
            if (spike_array(j) > 0 )
                plot([t(j) t(j)],[0 0.8*spike_array(j)],'k');
            end
        end
        axis([-before_stim/1000.0 (1000*response.stim_duration + after_stim)/1000.0 0 1]);
        axis off;
        hold off;
    end
end
psth = psth./ntrials;

subplot(4,1,3)
cla;
wind1 = hanning(31)/sum(hanning(31));   % 31 ms smoothing
smpsth = conv(psth,wind1);
plot(t,smpsth(16:length(smpsth)-15)*1000);
axis([-before_stim/1000.0 (1000*response.stim_duration + after_stim)/1000.0 0 1000*max(smpsth)]);
ylabel('Rate (spikes/s)');
xlabel('Time (s)');

% Add a nice title  
subplot(4,1,1);
titleName = sprintf('Stim %s %s(%s)', response.number, response.stim_class, response.stim_type);
if isfield(response, 'callid')
    titleName = strcat(titleName,sprintf(':%s', response.callid));
end
if isfield(response, 'birdid')
    titleName = strcat(titleName,sprintf('(%s)', response.birdid));
end
if isfield(response, 'callerAge') && isfield(response, 'stim_source_sex') && isfield(response, 'stim_source')
    titleName = strcat(titleName,sprintf('(%s %s %s)', response.callerAge, response.stim_source_sex, response.stim_source));
end
if isfield(response, 'distance')
    titleName = strcat(titleName,sprintf(' d=%s', response.distance));
end

titleName = strcat(titleName, sprintf(' z=%.1f p=%.3f', response.zscore, response.pvalue));

title(titleName);

%plot waveline and on
subplot(4,1,4);
wave=[zeros(before_stim/1000.0*samprate,1); sound_in; zeros(ceil(after_stim/1000*samprate),1)];
Max=0.1*ceil(10*max(wave));
plot(wave)
hold on
after_on = total_duration/1000*samprate - before_stim/1000*samprate - length(on);
on=[zeros(before_stim/1000.0*samprate,1); on; zeros(ceil(after_on),1)];
plot(Max*on, 'r');
hold on 
axis([0 length(on) -Max Max])
hold off;

