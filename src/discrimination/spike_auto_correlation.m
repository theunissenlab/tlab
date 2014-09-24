function outputs = spike_auto_correlation(birdname, brainregion, cellname, stimtype)


% Read input data from data base
[nfiles spike_times stim_len] = read_all_spikes(birdname, brainregion, cellname, stimtype);


% Initialize all return values to zero
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
auto_corr_fun = zeros(1,corr_time*10+1);

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
                        % break;    % delete this break to get auto_correlation. Keep it to get ISI
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
auto_corr_fun(1)=0;  % set to zero for displaying
subplot(1,1,1);
plot(x_time, auto_corr_fun);
%pause;

p1_100ms = sum(auto_corr_fun(1:10))./sum(auto_corr_fun);

outputs = struct('nfiles', nfiles, 'stimrate', stimrate, 'stdstimrate', stdstimrate, 'backrate', backrate, 'stdbackrate', stdbackrate, ...
    'avgzscore', avgzscore, 'stdzscore', stdzscore, 'avgdprime', avgdprime, 'stddprime', stddprime, ...
    'p1_100ms', p1_100ms );
    
return




