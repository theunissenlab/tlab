function [spikesfiltered, Spikeh, Hwidth] = gauss_filter_varying_window(spikes,  alpha_param, Fig)
%[spikesfiltered] = gauss_filter_varying_window(spikes, sd);
%this smooths with a gaussian window, with sd=sd to the furthest neighbor of psth (mean over all trials).

% spikesfiltered is the same size as spikes and is the gaussian filtered
% output
% Spikeh contains the distance to the nearest neighbor of each spike
% Hwidth is a structure that contains three vectors of the same length as 
% the number of spikes in the input spikes: Hwidth.trial is the trial id,
% Hwidth.timebin is the time bin position of each spike and Hwidth.hwidth
% is 1/2 the width of the gaussian each spike is convolved with.

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
    Spikeh = [];
    display('there are no spikes'); % the output is the unchanged input
    return;
end
if length(SpikeID)==1;
    spikesfiltered=ones(size(spikes))/size(spikes,2); % the input is a matrix (same size as input) of ones divided over the number of trials
    Spikeh = [];
    display('only one spike');
    return;
end

Spikeh = length(SpikeID);% This will contain the distance to the nearest neighbor of each spike
ss=0;
for nt=1:size(spikes,1)
    index=find(spikes(nt,:));% this is the bins in which there was at least a spike over all trials
    timevec = [0 index size(spikes,2)+1];% This is to deal with the first and last spike distances
    for a = 1:length(index) % iterate through the number of bins that contain spikes
        ss=ss+1;
        currenttime=index(a);% this is the time bin where is the current spike
        if spikes(nt,currenttime)>1% there is more than a spike in this time window
            Spikeh(ss) = 1;%This is the smallest distance we can handle? 1ms?
        else
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
end

ntrials=size(spikes,1);
temp=zeros(size(spikes));

if Fig
    F2=figure(2);
    for t=1:ntrials
        subplot(ntrials,1,t)
        cla
    end
end
Ymaxcc = nan(length(SpikeID),1);
Hwidth.trial = nan(length(SpikeID),1);
Hwidth.timebin = nan(length(SpikeID),1);
Hwidth.hwidth = nan(length(SpikeID),1);
for a = 1:length(SpikeID);
    hwidth=Spikeh(a)*alpha_param; % we want the number of points of the gaussian to be proportional to alpha so that the distance to the nearest spike is always = 1sd
    tempwin=gausswin(hwidth*2+1, alpha_param)/sum(gausswin(hwidth*2+1, alpha_param));
    tempgauss = tempwin';
    [currenttime,trial]=ind2sub(size(trans_spikes),SpikeID(a));
    Hwidth.trial(a) = trial;
    Hwidth.timebin(a) = currenttime;
    Hwidth.hwidth(a) = hwidth;
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
    temp(trial,startindex_temp:endindex_temp)=temp(trial,startindex_temp:endindex_temp) + spikes(trial,currenttime) * tempgauss(startindex_gauss:endindex_gauss);
    if Fig
        subplot(ntrials,1,trial)
        plot(startindex_temp:endindex_temp, tempgauss(startindex_gauss:endindex_gauss),'-r')
        Ymaxcc(a) = max(tempgauss(startindex_gauss:endindex_gauss));
        hold on
        plot([currenttime currenttime], [0 spikes(trial,currenttime)*0.1],'-k','LineWidth',1.5)
        hold on
        xlim([0 size(spikes,2)+1]);
        
    end
end
if Fig
    F2.Name='gaussian windows of all trials';
    subplot(ntrials,1,ntrials)
    xlabel('Time in ms')
%     Ymaxcc=nan(length(F2.Children),1);
%     for cc=1:length(F2.Children)
%         Ymaxcc(cc) = F2.Children(cc).YLim(2);
%     end
    YMAX=max(Ymaxcc);
    YMAX=YMAX+YMAX/5;
    if YMAX<0.2
        YMAX=0.2;
    end
    for cc=1:length(F2.Children)
        F2.Children(cc).YLim(2)=YMAX;
    end
    hold off
end
spikesfiltered=temp; % This is of the same size as the original spikes input
YMAX2 = max(max(temp));
YMAX2=YMAX2+YMAX2/5;

if Fig
    F3=figure(3);
    for t=1:ntrials
        subplot(ntrials,1,t)
        cla
        Spikes_local = find(spikes(t,:));
        for ss=1:length(Spikes_local)
            plot([Spikes_local(ss) Spikes_local(ss)],[0 YMAX2/5],'-k','LineWidth',1.5);
            hold on
        end
        plot((1:size(spikes,2)), spikesfiltered(t,:),'-r', 'LineWidth',1.5)
        xlim([0 size(spikes,2)+1])
        ylim([0 YMAX2])
    end
    F3.Name=sprintf('Gaussian filtered spike train Alpha=%d', alpha_param);
    subplot(ntrials,1,1)
    ylabel('spike rate spike/ms')
    subplot(ntrials,1,ntrials)
    xlabel('Time (ms)')
    hold off
    figure(4)
    cla
    [AX,H1,H2]=plotyy(1:size(spikes,2),mean(spikes,1),1:size(spikes,2),mean(spikesfiltered,1));
    ylabel(AX(1),'Average Spike Count')
    ylabel(AX(2), 'Filtered Spike rate')
    xlabel('Time in ms')
    xlim(AX(1),[0 size(spikes,2)+1])
    xlim(AX(2),[0 size(spikes,2)+1])
    ylim(AX(1),[0 YMAX2])
    ylim(AX(2),[0 YMAX2])
    set(H1, 'LineWidth', 1.5)
    set(H2, 'LineWidth', 1.5)
    set(AX(1), 'LineWidth', 1.5)
    set(AX(2), 'LineWidth', 1.5)
    title(sprintf('average result of the gaussian filtering over trials with alpha=%d',alpha_param));
end

end
