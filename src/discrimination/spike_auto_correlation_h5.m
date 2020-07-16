function outputs = spike_auto_correlation_h5(h5Path, stimType, perm)
% Calculates the inter-spike interval distribution.
% This is useful to check the quality of single units since one would
% expect to find deviations from exponential distribution because of
% refractory period. If perm is set to 'r' then only calculation, if set to
% 'a' then the code adds fields in the h5 to include outputs.

if nargin<3
    perm='r';
end

%% Read input data from data base
if strcmp(perm, 'r')
    unit = read_unit_h5file(h5Path, 'r');
else
    unit = read_unit_h5file(h5Path, 'a');
    h5 = h5utils();
    fid = h5.open(h5Path, 'a');
end

%% Select the protocol (SelX, Callx, Maskx, STRFx...)
%number of different protocols run for that unit, identify if the one
%you are asking for is existing in this unit and selecting this one in
%responses
nclasses = length(unit.classes);
classId = 0;
for ic=1:nclasses
    prot=unit.classes{ic};
    [path, Name, ext]=fileparts(unit.source_directory);
    if strcmp('WhiWhi4522M', Name) && strcmp('Site1', unit.site)
        if strcmp('Call1c', prot(1:6))
            classId = ic;
            break;
        end
    elseif strcmp(stimType, prot(1:4))
        classId = ic;
        break;
    end
end

if classId == 0
    fprintf(1, 'Warning: could not find stimType %s in h5 file %s\n', stimType, h5Path);
    fprintf(1, '\tAvailable options are:\n');
    for ic=1:nclasses
        fprintf(1,'\t\t%s\n', unit.classes{ic});
    end
    return;
end

responses = unit.class_responses.(unit.classes{ic});

%% Construct the input we need for all the following code
% This is the number of sound files played
nfiles = length(responses);

%This is the vector of the length of the different stims
stim_len = zeros (nfiles,1);
for nf=1:nfiles
    stim_len(nf) = 1000*responses{nf}.stim_duration;%here stim_len is in ms
end

% This is the cell array of the time of spike arrival
spike_times = cell(1, nfiles);
for nf = 1:nfiles
    Trials = responses{nf}.trials;
    nt = length(Trials);
    spike_times{nf} = cell(1, nt);
    for it = 1:nt
        trial = Trials{it};
        spike_times{nf}{it} = trial.spikeTimes.*1000.0;   % spike_times are in ms in neural_discrimination and in s in hdf5
    end
end


%% Initialize all return values to zero
stimrate = 0;
stdstimrate = 0;
backrate = 0;
stdbackrate = 0;
avgzscore = 0;
stdzscore = 0;
avgdprime = 0;
stddprime = 0;

tot_spikes = 0;
corr_time = 100;                        % Auto-correlation time in ms
auto_corr_fun = zeros(1,corr_time*10+1); % A bin size every 0.1 ms

if nfiles
    % Calculate firing rates, zscore and dprime
    [stimrate stdstimrate backrate stdbackrate avgzscore stdzscore avgdprime stddprime] = dprime_within(nfiles, spike_times, stim_len);

    % Calculate the spike time auto-correlation-function.
    for nfi=1:nfiles
        ntrials = length(spike_times{nfi});
        for it=1:ntrials
            spike_trials = spike_times{nfi}{it};
            ns = length(spike_trials);
            for is=1:ns
                if (spike_trials(is) < 0)
                    continue;
                end
                if (spike_trials(is) > stim_len(nfi))
                    break;
                end
                tot_spikes = tot_spikes+1;
                for is2=is+1:ns
                    if (spike_trials(is2) > stim_len(nfi))
                        break;
                    end
                    time_ind = round(10*(spike_trials(is2) - spike_trials(is))) + 1;
                    if time_ind <= corr_time*10+1
                        auto_corr_fun(time_ind) = auto_corr_fun(time_ind) + 1;
                        break;    % delete this break to get auto_correlation. Keep it to get ISI
                    else
                        break;
                    end
                end
            end
        end
    end
end
auto_corr_fun = auto_corr_fun ./ tot_spikes;

figure(1);
x_time = 0:0.1:corr_time;

%  auto_corr_fun(1)=0;  % if calculating a real auto-correlation function

subplot(1,1,1);
plot(x_time, auto_corr_fun);
xlabel('Time (ms)');
ylabel('Probability');
title('Inter Spike Interval Distribution');
%pause;

p1 = sum(auto_corr_fun(1:10))./sum(auto_corr_fun);
p1_theo = expcdf(0.001,1./mean(stimrate));

outputs = struct('nfiles', nfiles, 'stimrate', stimrate, 'stdstimrate', stdstimrate, 'backrate', backrate, 'stdbackrate', stdbackrate, ...
    'avgzscore', avgzscore, 'stdzscore', stdzscore, 'avgdprime', avgdprime, 'stddprime', stddprime, ...
    'p1', p1, 'p1_theo', p1_theo );

%% Write to file if requested
if strcmp(perm, 'a')
    pathPrefix = [sprintf('/extra_info/%s', prot) '/inter_spike_interval'];
    h5.set_ds(fid, pathPrefix, 'Spike_probability', auto_corr_fun);
    h5.set_ds(fid, pathPrefix, 'Interval_Inter_spike', x_time);
    h5.set_ds(fid, pathPrefix, 'Spike_proba_1ms', p1);
    h5.set_ds(fid, pathPrefix, 'Spike_proba_1ms_theo', p1_theo);
    h5.set_ds(fid, pathPrefix, 'Stim_rate', stimrate);
    h5.close(fid);
end
    
return




