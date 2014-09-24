function [talvecnewtotal, talvecoldtotal, rescaled_spikes]=time_rescaled_tals2(spiketrain, ratevec)
%[talvecnewtotal, talvecoldtotal]=time_rescaled_tals(spiketrain, ratevec);
%this gets time rescaled tals and original tals (tal=isi)
ntrials=size(spiketrain, 1);
difftimetotal=[];
meanratestotal=[];
talvecnewtotal=[];
talvecoldtotal=[];
burstvec=[.45 .3 .2];  
meanrate=mean(ratevec);
if meanrate==0
    display('no spikes in mean');
    rescaled_spikes=spiketrain;
    return
end
rescaled_spikes=zeros(size(spiketrain));
for a = 1:ntrials
    time_indices=find(spiketrain(a,:));
    extratimes=[];
    %the following takes care of bins with more than one spike
    for nspikes=2:4
        burstindex=find(spiketrain(a,:)==nspikes);
        for ns=1:nspikes-1
            extratimes=[extratimes burstindex+burstvec(nspikes-1)*ns];
        end
    end
    

    newtime=cumsum(ratevec)/meanrate;
    ntime_indices=newtime(time_indices);
    if ~isempty(extratimes)
        ntimeindicesextras=spline(1:length(ratevec), newtime, extratimes);
        ntime_indices=sort([ntime_indices ntimeindicesextras]);
    end
    talvecnew=diff(ntime_indices);
    talvecold=diff(sort([time_indices extratimes]));
    if ~isempty(ntime_indices)
        %this doesn't make a one if there is a zero.
        rescaled_spikes(a, round(ntime_indices(find(round(ntime_indices)))))=1;
    end
    talvecnewtotal=[talvecnewtotal talvecnew];
    talvecoldtotal=[talvecoldtotal talvecold];
end

