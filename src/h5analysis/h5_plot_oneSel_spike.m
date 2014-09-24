function h5_plot_one_spike(response,fsnip)
% plots the spectrogram, spike arrival times and average firing rate for
% response structure in our tdt hdf5 files. 
% fsnip is the pointer to the snipet file or -1 if no snipets should be
% plotted.

before_stim = 2000;                     % Time before stimulus
total_duration = 6000;                      % Time after stimulus
nsnip = 18;                              % Length of a snippet - the length is set in RPVDs macro used for TDT.
% This value of 20 was taken from 32Channel2EC.rcx in
% Masker_Multielectrode_FJN.  THe snipped length is set in the
% Spike_Store_MC macro.


% Read the stim wave files on the cluster on a local mac machine.
stim_name = response.tdt_wavfile;

if ismac()
    [status username] = system('who am i');
    if strcmp(strtok(username), 'frederictheunissen')
        if strncmp('/auto/fdata/solveig',stim_name, 19)
            stim_name = strcat('/Users/frederictheunissen/Documents/Data/solveig', stim_name(20:end));
        elseif strncmp('/auto/fdata/julie',stim_name, 17)
            stim_name = strcat('/Users/frederictheunissen/Documents/Data/Julie', stim_name(18:end));
        end
    elseif strcmp(strtok(username), 'elie')
        if strncmp('/auto/fdata/solveig',stim_name, 19)
            stim_name = strcat('/Users/frederictheunissen/Documents/Data/solveig', stim_name(20:end));
        elseif strncmp('/auto/fdata/julie',stim_name, 17)
            stim_name = strcat('/Users/elie/Documents/MATLAB/data', stim_name(18:end));
        end
    elseif strcmp(strtok(username), 'Solveig')
        if strncmp('/auto/fdata/solveig',stim_name, 19)
            stim_name = strcat('/Users/Solveig/PhD', stim_name(20:end));
        end
    end
end
[sound_in, samprate] = wavread(stim_name);

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
%figure(); % this is changed from h5_plot_one_spike to avoid using the same figure 1 for all plots
h = gcf();
set(h,'PaperPosition',[0.25 2.5 3.75 8]); % This is for hard copy
set(h,'Position', [50 50 1000 500]);

% Calculate and plot the spectrogram   
subplot(3,1,1); 
cla;
[s, to, fo, pg] = GaussianSpectrum(sound_in, increment, winLength, samprate); 
logB = 20*log10(abs(s));
maxB = max(max(logB));
minB = maxB-DBNOISE;            

imagesc(to,fo,logB);          % to is in seconds
axis xy;
caxis('manual');
caxis([minB maxB]); 

after_stim = total_duration - before_stim - response.stim_duration;
v_axis(1) = -before_stim/1000.0;
v_axis(2) = (response.stim_duration + after_stim)/1000.0;
v_axis(3)=f_low; 
v_axis(4)=f_high;

axis(v_axis);                                
xlabel('time (ms)'), ylabel('Frequency');
cmap = spec_cmap();
colormap(cmap);

% plot spike array and calculate average
ntimebins = before_stim + response.stim_duration + after_stim; % Number of time bins in ms
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
        subplot('position',[0.13 0.66-it*(0.3/nplot) 0.775 (1-0.05)*0.3/nplot]);
        cla;
        hold on;
        for j=1:ntimebins
            if (spike_array(j) > 0 )
                plot([t(j) t(j)],[0 spike_array(j)],'k');
            end
        end
        axis([-before_stim/1000.0 (response.stim_duration + after_stim)/1000.0 0 1]);
        axis off;
        hold off;
    end
end
psth = psth./ntrials;

subplot(3,1,3)
cla;
wind1 = hanning(31)/sum(hanning(31));   % 31 ms smoothing
smpsth = conv(psth,wind1);
plot(t,smpsth(16:length(smpsth)-15)*1000);
axis([-before_stim/1000.0 (response.stim_duration + after_stim)/1000.0 0 1000*max(smpsth)]);
ylabel('Rate (spikes/s)');
xlabel('Time (s)');

% Add a nice title  
subplot(3,1,1);
titleName = sprintf('Stim %s %s(%s)', response.number, response.stim_class, response.stim_type);
if isfield(response, 'callid')
    titleName = strcat(titleName,sprintf(':%s', response.callid));
end
if isfield(response, 'birdid')
    titleName = strcat(titleName,sprintf('(%s)', response.birdid));
end
if isfield(response, 'callnum')
    titleName = strcat(titleName,sprintf('(call%s)', response.callnum));
end
if isfield(response, 'callerAge') && isfield(response, 'stim_source_sex') && isfield(response, 'stim_source')
    titleName = strcat(titleName,sprintf('(%s %s %s)', response.callerAge, response.stim_source_sex, response.stim_source));
end
if isfield(response, 'distance')
    titleName = strcat(titleName,sprintf(' d=%s', response.distance));
end

titleName = strcat(titleName, sprintf(' z=%.1f p=%.3f', response.zscore, response.pvalue));

title(titleName);

if fsnip ~= -1
    figure(2);
    for it=1:ntrials
        trial = response.trials{it};
        spike_id_trials = trial.spikeIds;
        ns = length(spike_id_trials);
        for is=1:ns
            fseek(fsnip, (spike_id_trials(is)-1)*4*nsnip, 'bof');
            snip = fread(fsnip, nsnip, 'float32');
            plot(snip,'k');
            hold on;
        end
    end
    hold off
end
