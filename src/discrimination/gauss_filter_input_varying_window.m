function [spikesfiltered, Hwidth_out] = gauss_filter_input_varying_window(spikes, Hwidth_in, Fig)
%this smooths with a gaussian window, using the width that is input with Hwidth.

% spikesfiltered is the same size as spikes and is the gaussian filtered
% output
% Hwidth is a structure that contains three vectors of the same length:
% Hwidth.trial is the trial id,
% Hwidth.timebin is the time bin position of each spike width previously
% calculated and Hwidth.hwidth is 1/2 the width of the gaussian each spike
% should be convolved with.
% Hwidth.alpha_param is the alpha parameter that was previously used to calculate the
% width of the gaussians

% spikes is an n*t matrix with n the number of responses obtained to the...
...same stimulus presentation (can be one), and t the number of time bins


if nargin<3
    Fig=1; %set to 1 to see figures
end
trans_spikes=spikes';
SpikeID =  find(trans_spikes);% This is the linear index of the spikes in the input matrix spikes

if isempty(SpikeID)
    spikesfiltered=spikes;
    Hwidth_out = [];
    fprintf(1,'there are no spikes\n'); % the output is the unchanged input
    return;
end
if length(SpikeID)==1
    spikesfiltered=ones(size(spikes))/size(spikes,2); % the input is a matrix (same size as input) of ones divided over the number of trials
    Hwidth_out = [];
    fprintf(1,'only one spike\n');
    return;
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
Hwidth_out.trial = nan(length(SpikeID),1);
Hwidth_out.timebin = nan(length(SpikeID),1);
Hwidth_out.hwidth = nan(length(SpikeID),1);
for SI = 1:length(SpikeID)
    [currenttime,trial]=ind2sub(size(trans_spikes),SpikeID(SI));
    Ind4Trial = find(Hwidth_in.trial == trial);
    TimeBin4Trial = Hwidth_in.timebin(Ind4Trial);
    Hwidth4Trial = Hwidth_in.hwidth(Ind4Trial);
    
    if sum(TimeBin4Trial==currenttime)
        hwidth = Hwidth4Trial(find(TimeBin4Trial==currenttime));
    else
        Latertime = find(TimeBin4Trial>currenttime);
        Beforetime = find(TimeBin4Trial<currenttime);
        if isempty(Latertime) 
            hwidth = Hwidth4Trial(Beforetime(1));
        elseif isempty(Beforetime)
            hwidth = Hwidth4Trial(Latertime(1));
        else
            y1 = Hwidth4Trial(Beforetime(1));
            y2 = Hwidth4Trial(Latertime(1));
            x1= TimeBin4Trial(Beforetime(1));
            x2 = TimeBin4Trial(Latertime(1));
            a = (y1-y2)/(x1-x2);
            b = y1-a*x1;
            hwidth = round(a*currenttime + b);
        end
    end
    
    tempwin=gausswin(hwidth*2+1, Hwidth_in.alpha_param)/sum(gausswin(hwidth*2+1, Hwidth_in.alpha_param));
    tempgauss = tempwin';
    
    Hwidth_out.timebin(SI) = currenttime;
    Hwidth_out.hwidth(SI) = hwidth;
    if currenttime-(hwidth+1)<=0 % The gaussian left half width is larger than the distance of the spike to the begining of the spike train
        startindex_temp = 1;
        startindex_gauss = hwidth+1 -currenttime+1;
    else
        startindex_gauss=1;
        startindex_temp = currenttime-hwidth;
    end
    if currenttime+hwidth>=size(spikes,2) % The gaussian right half width is larger than the distance of the spike to the end of the spike train
        endindex_temp=size(spikes,2);
        endindex_gauss = hwidth+1 + size(spikes,2) -currenttime;
    else
        endindex_temp = currenttime + hwidth;
        endindex_gauss = size(tempgauss,2);
    end
    temp(trial,startindex_temp:endindex_temp)=temp(trial,startindex_temp:endindex_temp) + spikes(trial,currenttime) * tempgauss(startindex_gauss:endindex_gauss);
    if Fig
        subplot(ntrials,1,trial)
        plot(startindex_temp:endindex_temp, spikes(trial,currenttime) *tempgauss(startindex_gauss:endindex_gauss),'-r')
        Ymaxcc(SI) = max(spikes(trial,currenttime) *tempgauss(startindex_gauss:endindex_gauss));
        hold on
        plot([currenttime currenttime], [0 spikes(trial,currenttime)*0.1],'-k','LineWidth',1.5)
        hold on
        xlim([0 size(spikes,2)+1]);
        
    end
end
if Fig
    F2.Name='gaussian windows of all trials';
    subplot(ntrials,1,ntrials)
    xlabel('Time in time bins')
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
    F3.Name=sprintf('Gaussian filtered spike train Alpha=%d', Hwidth_in.alpha_param);
    subplot(ntrials,1,1)
    ylabel('spike rate spike/timebin')
    subplot(ntrials,1,ntrials)
    xlabel('Time (time bins)')
    hold off
    figure(4)
    cla
    [AX,H1,H2]=plotyy(1:size(spikes,2),mean(spikes,1),1:size(spikes,2),mean(spikesfiltered,1));
    ylabel(AX(1),'Average Spike pattern')
    ylabel(AX(2), 'Filtered Spike pattern')
    xlabel('Time in ms')
    xlim(AX(1),[0 size(spikes,2)+1])
    xlim(AX(2),[0 size(spikes,2)+1])
    ylim(AX(1),[0 YMAX2])
    ylim(AX(2),[0 YMAX2])
    set(H1, 'LineWidth', 1.5)
    set(H2, 'LineWidth', 1.5)
    set(AX(1), 'LineWidth', 1.5)
    set(AX(2), 'LineWidth', 1.5)
    title(sprintf('average result of the gaussian filtering over trials with alpha=%d',Hwidth_in.alpha_param));
end

end
