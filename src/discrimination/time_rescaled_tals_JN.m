function [talvecnewtotal, talvecoldtotal, rescaled_spikes, smoothedspike]=time_rescaled_tals_JN(spiketrain, sd);
%time_rescaled_tals_JN(spiketrain);
%this performs time_rescaled_tals2 on each trial of spike train using mean rate as the gauss_filter_psth_near
%on the mean of all trials but that one.
ntrials=size(spiketrain,1);
talvecnewtotal=[];
talvecoldtotal=[];
rescaled_spikes=zeros(size(spiketrain));
spikefilt=gauss_filter_vwbytrial(spiketrain, sd);
%smoothedspike=gauss_filter_vwbytrial(spiketrain, 2.5);
for nt=1:ntrials
    meanrate=mean(spikefilt([1:nt-1 nt+1:end],:));
    %meanrate=mean(tempspike);

    [talvecnew, talvecold, rescaled_spiketrain]=time_rescaled_tals2(spiketrain(nt,:), meanrate);
    rescaled_spikes(nt,:)=rescaled_spiketrain;
    talvecnewtotal=[talvecnewtotal talvecnew];
    talvecoldtotal=[talvecoldtotal talvecold];
end