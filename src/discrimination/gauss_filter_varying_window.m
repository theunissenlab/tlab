function [spikesfiltered, Spikeh] = gauss_filter_varying_window(spikes,  alpha_param, Fig)
%[spikesfiltered] = gauss_filter_varying_window(spikes, sd);
%this smooths with a gaussian window, with sd=sd to the furthest neighbor of psth (mean over all trials).

% spikes is an n*t matrix with n the number of responses obtained to the...
...same stimulus presentation (can be one), and t the number of time bins in ms

%timevec is a row vector of times of all spikes present per trial
%and includes 0 and size(spikes,2) at beginning and end of each trial
%to identify first and last spike

if nargin<2
    alpha_param=3;
end
if nargin<3
    Fig=1; %set to 1 to see figures
end
trans_spikes=spikes';
SpikeID =  find(trans_spikes);% This is the linear index of the spikes in the input matrix spikes

if isempty(SpikeID);
    spikesfiltered=spikes;
    display('there are no spikes'); % the output is the unchanged input
    return;
end
if length(SpikeID)==1;
    spikesfiltered=ones(size(spikes))/size(spikes,2); % the input is a matrix (same size as input) of ones divided over the number of trials
    display('only one spike');
    return;
end

Spikeh = length(SpikeID);% This will the contain the distance to the nearest neighbor of each spike
ss=0;
for nt=1:size(spikes,1)
    index=find(spikes(nt,:));% this is the bins in which there was at least a spike over all trials
    timevec = [0 index size(spikes,2)+1];% This is to deal with the first and last spike distances
    for a = 1:length(index) % iterate through the number of bins that contain spikes
        ss=ss+1;
        currenttime=index(a);% this is the time bin where is the current spike
        timeindex=find(timevec==index(a));% This is the position of the spike in timevec
        futuretime=timevec(timeindex+1);% this is the arrival time of the following spike
        pasttime=timevec(timeindex-1);% this is the arrival time of the preceding spike
        after=futuretime-currenttime;% This is the time between the spike and the future spike
        before=currenttime-pasttime;% This is the time between the spike and the preceeding spike
        if futuretime==size(spikes,2)+1 % The current spike was the last spike
            after=after*2+1; % More than double the time after the spike
        elseif pasttime==0 % The current spike was the first spike of the trial
            before=before*2+1; % more than double the time before the spike
        end
        Spikeh(ss)=min([before after]);% This is the distance in number of bins to the nearest neighbor spike
    end
end

ntrials=size(spikes,1);
temp=zeros(size(spikes));

if Fig
    figure(2)
    for t=1:ntrials
        subplot(ntrials,1,t)
        cla
    end
end
for a = 1:length(SpikeID);
    hwidth=Spikeh(a)*alpha_param; % we want the number of points of the gaussian to be proportional to alpha so that the distance to the nearest spike is always = 1sd
    tempwin=gausswin(hwidth*2+1, alpha_param)/sum(gausswin(hwidth*2+1, alpha_param));
    tempgauss = tempwin';
    [currenttime,trial]=ind2sub(size(trans_spikes),SpikeID(a));
    if currenttime-(hwidth+1)<=0 % The gaussian left half width is larger than the distance of the spike to the begining of the spike train
        startindex_temp = 1;
        startindex_gauss = hwidth+1 -currenttime+1;
    else
        startindex_gauss=1;
        startindex_temp = currenttime-hwidth;
    end
    if currenttime+hwidth>=size(spikes,2);% The gaussian right half width is larger than the distance of the spike to the end of the spike train
        endindex_temp=size(spikes,2);
        endindex_gauss = hwidth+1 + size(spikes,2) -currenttime;
    else
        endindex_temp = currenttime + hwidth;
        endindex_gauss = size(tempgauss,2);
    end
    temp(trial,startindex_temp:endindex_temp)=temp(trial,startindex_temp:endindex_temp)+tempgauss(startindex_gauss:endindex_gauss);
    if Fig
        subplot(ntrials,1,trial)
        plot(startindex_temp:endindex_temp, tempgauss(startindex_gauss:endindex_gauss),'-r')
        hold on
        plot([currenttime currenttime], [0 spikes(trial,currenttime)*0.1],'-k')
        hold on
        xlim([0 size(spikes,2)+1]);
        
    end
end
if Fig
    figure(2)
    xlabel('Time in ms')
    title('gaussian windows of all trials')
    hold off
end
spikesfiltered=temp; % This is of the same size as the original spikes input

if Fig
    figure(3)
    for t=1:ntrials
        subplot(ntrials,1,t)
        cla
        Spikes_local = find(spikes(t,:));
        for ss=1:length(Spikes_local)
            plot([Spikes_local(ss) Spikes_local(ss)],[0 0.25],'-k');
            hold on
        end
        plot((1:size(spikes,2)), spikesfiltered(t,:),'-r')
        title(sprintf('Gaussian filtered spike train %d, Alpha=%d',t, alpha_param));
        xlabel('Time in ms')
        ylabel('spike rate spike/ms')
        xlim([0 size(spikes,2)+1])
    end
    hold off
    figure(4)
    cla
    [AX,~]=plotyy(1:size(spikes,2),mean(spikes,1),1:size(spikes,2),mean(spikesfiltered,1));
    ylabel(AX(1),'Average Spike Count')
    ylabel(AX(2), 'Filtered Spike rate')
    xlabel('Time in ms')
    xlim(AX(1),[0 size(spikes,2)+1])
    xlim(AX(2),[0 size(spikes,2)+1])
    title(sprintf('average result of the gaussian filtering over trials with alpha=%d',alpha_param));
end

end
