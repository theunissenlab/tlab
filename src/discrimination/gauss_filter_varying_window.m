function [spikesfiltered, h, index] = gauss_filter_varying_window(spikes,  sd);
%[spikesfiltered] = gauss_filter_varying_window(spikes, sd);
%this smooths with a gaussian window, with sd=sd to the furthest neighbor of psth (mean over all trials).

%timevec is a row vector of times of all spikes present in all trials
%and includes 0 and length(spiketrain) at beginning and end of each trial
%all trials are in a big long row.

spiketrain=mean(spikes,1);
timevec=[];
for a=1:size(spikes,1)
timevec=[timevec 0 find(spikes(a,:)) size(spikes,2)+1];
end
index=find(spiketrain);
if isempty(index);
    spikesfiltered=spikes;
    display('there are no spikes');
    return;
end
if length(index)==1;
    spikesfiltered=ones(size(spikes))/size(spikes,2);
    display('only one spike');
    return;
end
    h=zeros(1,length(index));
for a = 1:length(index)
    timeindex=find(timevec==index(a));
    currenttime=timevec(timeindex);
    futuretime=timevec(timeindex+1);
    pasttime=timevec(timeindex-1);
    after=futuretime-currenttime;
    before=currenttime-pasttime;
    if ~isempty(find(futuretime==size(spikes,2)+1))
        after=after*2+1;
    elseif ~isempty(find(pasttime==0))
        before=before*2+1;
    end
    h(a)=max([before after]);
    

end
h=ceil(h/2)+1;



ntrials=size(spikes,1);
temp=zeros(size(spikes));

%winvec=zeros(1,length(spiketrain));
for a = 1:length(index);
    hwidth=h(a);
    tempwin=gausswin(hwidth*2+1, sd)/sum(gausswin(hwidth*2+1, sd));
    tempgauss=repmat(tempwin',ntrials,1).*repmat(spikes(:,index(a)),1, length(tempwin));
    if index(a)-hwidth<=0
        endindex=round(min([index(a)+hwidth length(temp)]));
        templength=endindex;
        temp(:,1:endindex)=temp(:,1:endindex)+tempgauss(:,end-templength+1:end);
        %winvec(1:index(a)+hwidth)=winvec(1:index(a)+hwidth)+tempwin(hwidth-index(a)+2:end);
    elseif index(a)+hwidth>length(spiketrain);
        startindex=round(max([index(a)-hwidth 1]));
        templength=size(temp,2)-startindex+1;
        temp(:,startindex:end)=temp(:,startindex:end)+tempgauss(:,1:templength);
        % winvec(index(a)-hwidth:end)=winvec(index(a)-hwidth:end)+tempwin(1:end-(index(a)+hwidth-length(spiketrain)));

    else
        endindex=round(min([index(a)+hwidth length(temp)]));
        startindex=round(max([index(a)-hwidth 1]));
        temp(:,startindex:endindex)=temp(:,startindex:endindex)+tempgauss;
        % winvec(index(a)-hwidth:index(a)+hwidth)=winvec(index(a)-hwidth:index(a)+hwidth)+tempwin;
    end
end

spikesfiltered=temp;
