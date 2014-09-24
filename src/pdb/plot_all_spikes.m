function plot_all_spikes(birdname, brainregion, cellname, stimtype)
% Plots spectrogram and spike rasters for all stimuli of one type recorded
% at one site.

[nfiles, spike_times, stimlen, stimname] = read_all_spikes(birdname, brainregion, cellname, stimtype);

for nfi=1:nfiles
    hf = plot_one_spike(spike_times{nfi}, stimlen(nfi), stimname{nfi});
    subplot(3,1,1);
    title(sprintf('%s %s %s %s %s', birdname, brainregion, cellname, stimtype, stimname{nfi})); 
    %pause();
end

return




