function spikesfiltered=gauss_filter_vwbytrial(spiketrain, sd);
spikesfiltered=zeros(size(spiketrain));
for ntrials=1:size(spiketrain,1)
    spikesfiltered(ntrials,:)=gauss_filter_varying_window(spiketrain(ntrials,:), sd);
end

spikesfiltered(find(spikesfiltered<2e-6))=2e-6;