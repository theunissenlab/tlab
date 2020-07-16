function spiketrain = spike_times_to_train(nfiles, spike_times, stim_len, onset_time)
% Converts spike arrival times given in spike_times to an time-series array of zeros
% and ones.  Spike times are assumed to be in ms and the time-series will
% be sampled at 1kHz.  The time-series from the different files are also
% concatenated into a single long time-series.  Only the spike times during
% the stimulus are preserved and the initial spikes in onset_time are
% ignored.  Set onset_time to zero to include all spikes.

spiketrain=[];

for nf=1:nfiles
    nt = length(spike_times{nf});
    spike = zeros(nt, round(stim_len(nf) - onset_time));

    for it=1:nt
        time_ind = round(spike_times{nf}{it}-onset_time);
        nspikes = length(time_ind);
        for is=1:nspikes
            if ( (time_ind(is) > 0) && (time_ind(is) <= round(stim_len(nf) - onset_time)) )
                spike(it,time_ind(is)) = 1;
            end
        end
    end
    if nf > 1
        if size(spike,1) ~=size(spiketrain,1);
            display(['spike' num2str(nf) ' skipped']);
        else
            spiketrain=[spiketrain spike];
        end
    else
        spiketrain=[spiketrain spike];
    end
end

return;