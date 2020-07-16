function [Spikesfilt, Hwidth] = gauss_filter_varying_window(Spikes, Alpha_param, Kth_neigh, Fig)

% Function to smooth spike rasters with a gaussian window, with sd=sd to the kth neighbor, in order to obtain estimate of
% time varying firing rate.

% Requires:
% Spikes        is an n*t matrix with n the number of responses (trials)
%               obtained to the same stimulus presentation (can be one),
%               and t the number of time bins in ms. Spikes can also be the
%               PSTH or concatenated spike pattern obtained from several
%               presentations of the same stimulus.


% Optional Input parameters
% Alpha_param   Alpha_param corresponds to the Alpha parameter in the
%               gausswin function of matlab. It is defined as the
%               reciprocal of the standard deviation and is a measure of
%               the width of its Fourier Transform. As alpha_param
%               increases, the width of the window will decrease.
%               If omitted, alpha_param is 3. 

% Kth_neigh     is the neighbor spike that determine the width of the
%               gaussian window. Kth_neigh is vector of n floating points
%               (one value for each row in the matrix spikes). The
%               default neighbor is set to the first nearest neighbor (1).
%               Note that the best value is 1/2 number of trials for
%               concatenated spike patterns of ntrials.

% Fig           a logical flag to see intermediare and output figures.
%               Set to 1 by default.


% Returns:
% Spikesfilt    is the same size as spikes and is the gaussian filtered
%               output of Spikes


% Hwidth        a structure that contains three vectors of length the number
%               of bins containing at least one spike in the input matrix Spikes.
%                       Hwidth.trial is the row/trial ID of the bin
%                       Hwidth.timebin is the time bin position of each
%                       bin containing at least a spike
%                       Hwidth.hwidth is 1/2 the width of the gaussian each
%                       spike in that bin is convolved with.

%% Set default input values
if nargin<2
    Alpha_param=3;
end
if nargin<3
    Kth_neigh = ones(size(Spikes,1),1);% The first neighbor spike is set...
    ... as default.
end
if nargin<4
    Fig=1; %set to 1 to see figures
end

%% Sanitary check of the input
% Calculate the linear index of the bins containing sipke(s) in the input matrix Spikes
% This indexation is used to easily refer to bins containing spikes.
Trans_spikes=Spikes';
SpikedBinID =  find(Trans_spikes);

if isempty(SpikedBinID)
    Spikesfilt=Spikes;
    Dist2KthSpike = [];
    HalfWidth = [];
    fprintf(1,'there are no spikes\n'); % the output is the unchanged input
    return;
end
if length(SpikedBinID)==1
    Spikesfilt=ones(size(Spikes))/size(Spikes,2); % the output is a matrix (same size as input) of ones divided over the number of rows/trials
    Dist2KthSpike = [];
    HalfWidth = [];
    fprintf(1,'only one spike\n');
    return;
end

%% Distance to the Kth neighbor spike
% This first part of the code identify the shortest distance in number of
% bins to the Kth nearest neighbor spike for all the spikes in each row of
% the matrix Spikes. Each row is treated independantly of the others.
% Spikes contained in the same bin have the same distance to the
% Kth nearest neighbor spike. This distance is reported in Dist2KthSpike that
% follows the linear indexation of bins containing spikes.

% Initializing the output vector
Dist2KthSpike = nan(1,length(SpikedBinID));

% Looping through rows/trials in Spikes and through the bins that contain spikes 
ss=0; % This index keeps track of the indices of the bins containing spikes
for nt=1:size(Spikes,1)
    Index=find(Spikes(nt,:));% this vector contains the indices of the bins in which there was at least a spike in this row
    for a = 1:length(Index)
        ss=ss+1;
        
        % The strategy here is to find the bin containing the Kth nearest
        % neighbor spike between past and future bins. Because Kth can be a
        % floating number, we are calculating the distance to both
        % the floor(Kth_neigh) and ceil(Kth_neigh) neighbor spikes and then obtain the
        % exact distance to the Kth neighbor by perfoming a linear
        % interpolation.
        
        % First find the bins containing the Kth future spikes
        AfterSpikes = cumsum(Spikes(nt,Index(a):end))-1; % cumulative number of other spikes in the present bin and in the following bins 
        AfterSpikesCeil = find(AfterSpikes >= (ceil(Kth_neigh(nt))), 1)-1; % Distance (in # of bins) to the bin in the future that contains the ceil(Kth_neigh) nearest neighbor spike 
        AfterSpikesFloor = find(AfterSpikes >= (floor(Kth_neigh(nt))), 1)-1; % Distance (in # of bins) to the bin in the future that contains the floor(Kth_neigh) nearest neighbor spike. ceil(Kth) and floor(Kth) are different if Kth is not a whole number
        % Second find the bins containing the Kth past spikes
        BeforeSpikes = cumsum(Spikes(nt,flip(1:Index(a))))-1; % cumulative number of other spikes in the present bin and in the precedent bins 
        BeforeSpikesCeil = find(BeforeSpikes >= (ceil(Kth_neigh(nt))), 1)-1; % Distance (in # of bins) to the bin in the past that contains the ceil(Kth_neigh) nearest neighbor spike 
        BeforeSpikesFloor = find(BeforeSpikes >= (floor(Kth_neigh(nt))), 1)-1; % Distance (in # of bins) to the bin in the past that contains the floor(Kth_neigh) nearest neighbor spike 
        
        % Treat the case where there are no Kth_neigh spikes in the future and
        % calculate the distance to the Kth_neigh future spike
        if isempty(AfterSpikesCeil) % There is no ceil(Kth_neigh) spikes in the future (not enough spikes after the current one)
            if isempty(AfterSpikesFloor) % And there is not even floor(Kth_neigh) spikes in the future (not enough spikes after the current one)
                After = [];
            else % There is at least floor(Kth_neigh) spikes in the future
                After = AfterSpikesFloor; % Consider the distance to this floor(Kth_neigh) spike as our best estimate of the distance to the Kth_neigh spike in the future
            end
        elseif Kth_neigh(nt)==fix(Kth_neigh(nt)) % Kth_neigh is a whole number, afterSpikesCeil = afterSpikesFloor 
            After = AfterSpikesFloor;
        else % Linear interpolation between bin positions of floor(Kth_neigh) and ceil(Kth_neigh) nearest neighbor spike if Kth_neigh is not a whole number
            After = AfterSpikesFloor + (AfterSpikesCeil - AfterSpikesFloor) * (Kth_neigh(nt) - floor(Kth_neigh(nt))) / (ceil(Kth_neigh(nt)) - floor(Kth_neigh(nt))) ;
        end
        
        % Treat the case where there are no Kth_neigh spikes in the past and
        % calculate the distance to the Kth_neigh past spike
        if isempty(BeforeSpikesCeil) % There is no ceil(Kth_neigh) spikes in the past (not enough spikes before the current one)
            if isempty(BeforeSpikesFloor) % There is not even floor(Kth_neigh) spikes in the past (not enough spikes before the current one)
                Before = [];
            else % There is at least floor(Kth_neigh) spikes in the past
                Before = BeforeSpikesFloor; % Consider the distance to this floor(Kth_neigh) spike as our best estimate of the distance to the Kth_neigh spike in the past
            end
        elseif Kth_neigh(nt)==fix(Kth_neigh(nt)) % Kth_neigh is a whole number, beforeSpikesCeil = beforeSpikesFloor 
            Before = BeforeSpikesFloor;
        else % Linear interpolation between bin positions of floor(Kth_neigh) and ceil(Kth_neigh) nearest neighbor spike if Kth_neigh is not a whole number
            Before = BeforeSpikesFloor + (BeforeSpikesCeil - BeforeSpikesFloor) * (Kth_neigh(nt) - floor(Kth_neigh(nt))) / (ceil(Kth_neigh(nt)) - floor(Kth_neigh(nt))) ;
        end
        
        % Treat the case when there is no Kth_neigh spikes both in the
        % future and in the past.Then set the distance to the distance to
        % the end or begining of the trial/row, whichever is the smallest
        if isempty(Before) && isempty(After)
            Dist2KthSpike(ss) = min([Index(a) size(Spikes,2) - Index(a)]);
        else
            % Find the minimum distance in number of bins to the Kth neighbor spike
            Dist2KthSpike(ss) = min([Before After]);
        end
        
        % Treat the case when several spikes are in the same bin and
        % Kth_neigh = 0. We want the minimum distance to be 0.5 and not 0
        % so the spikes are still convolved with a gaussian
        if Dist2KthSpike(ss) == 0
            Dist2KthSpike(ss) = 0.5;
        end
        
    end
end


%% Convolve spikes with gaussian windows of width proportional to the distance to the Kth_neigh nearest spike
% Extract the number of rows/trials in Spikes
Ntrials=size(Spikes,1);


% Initialize a figure that will plot the gaussian windows used for each
% row/trial
if Fig
    Ymaxcc = nan(length(SpikedBinID),1); % This will be used to scale the figure to the max of the gaussians
    F2=figure(2);
    for t=1:Ntrials
        subplot(Ntrials,1,t)
        cla
    end
end

% Initialize output variables
Spikesfilt = zeros(size(Spikes)); % This matrix will contain the gaussian filtered spike patterns
Hwidth.trial = nan(length(SpikedBinID),1); % This structure contains the information to identify the width of the gaussian with which spikes are convolved
Hwidth.timebin = nan(length(SpikedBinID),1);
Hwidth.halfwidth = nan(length(SpikedBinID),1);


% Now loop through Bins containing spike(s) and do the convolution
for a = 1:length(SpikedBinID)
    % set the width of the gaussian
    HalfWidth=round(Dist2KthSpike(a)*Alpha_param); % we want the number of points of the gaussian to be proportional to alpha so that the distance to the nearest spike is always = 1sd
    
    % calculate the gaussian window
    Tempwin=gausswin(HalfWidth*2+1, Alpha_param)/sum(gausswin(HalfWidth*2+1, Alpha_param));
    Tempgauss = Tempwin';
    
    % retrieve the bin position using its linear index and store position
    % and width of the gaussian
    [Currenttime,Trial]=ind2sub(size(Trans_spikes),SpikedBinID(a));
    Hwidth.trial(a) = Trial;
    Hwidth.timebin(a) = Currenttime;
    Hwidth.halfwidth(a) = HalfWidth;
    
    % Deal with cases where the width of the gaussian goes beyong the spike
    % train boundaries and set the begining and end indices where the
    % gaussian should be added to final gaussian filtered spike train(s)
    if Currenttime-(HalfWidth+1)<=0 % The gaussian left half width is larger than the distance of the spike to the begining of the spike train
        Startindex_temp = 1;
        Startindex_gauss = HalfWidth+1 -Currenttime+1;
    else
        Startindex_gauss=1;
        Startindex_temp = Currenttime-HalfWidth;
    end
    if Currenttime+HalfWidth>=size(Spikes,2) % The gaussian right half width is larger than the distance of the spike to the end of the spike train
        Endindex_temp=size(Spikes,2);
        Endindex_gauss = HalfWidth+1 + size(Spikes,2) -Currenttime;
    else
        Endindex_temp = Currenttime + HalfWidth;
        Endindex_gauss = size(Tempgauss,2);
    end
    
    % Convolve the gaussian window with the spike(s) and add to the output
    % matrix
    Spikesfilt(Trial,Startindex_temp:Endindex_temp)=Spikesfilt(Trial,Startindex_temp:Endindex_temp) + Spikes(Trial,Currenttime) * Tempgauss(Startindex_gauss:Endindex_gauss);
    
    % Fill in figure 2 with the gaussian window used for that bin
    if Fig
        subplot(Ntrials,1,Trial)
        plot(Startindex_temp:Endindex_temp, Spikes(Trial,Currenttime) *Tempgauss(Startindex_gauss:Endindex_gauss),'-r')
        Ymaxcc(a) = max(Spikes(Trial,Currenttime) *Tempgauss(Startindex_gauss:Endindex_gauss));
        hold on
        plot([Currenttime Currenttime], [0 Spikes(Trial,Currenttime)*0.1],'-k','LineWidth',1.5)
        hold on
        xlim([0 size(Spikes,2)+1]);
        
    end
end


% Fine tune figure 2 (legend....)
if Fig
    F2.Name='gaussian windows of all trials';
    subplot(Ntrials,1,Ntrials)
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


%% Plot a figure that gives the gaussian filtered spike pattern obtained for
% each row/trial in Spikes
if Fig
    YMAX2 = max(max(Spikesfilt));
    YMAX2=YMAX2+YMAX2/5;
    F3=figure(3);
    for t=1:Ntrials
        subplot(Ntrials,1,t)
        cla
        Spikes_local = find(Spikes(t,:));
        for ss=1:length(Spikes_local)
            plot([Spikes_local(ss) Spikes_local(ss)],[0 YMAX2/5],'-k','LineWidth',1.5);
            hold on
        end
        plot((1:size(Spikes,2)), Spikesfilt(t,:),'-r', 'LineWidth',1.5)
        xlim([0 size(Spikes,2)+1])
        ylim([0 YMAX2])
    end
    F3.Name=sprintf('Gaussian filtered spike train Alpha=%d', Alpha_param);
    subplot(Ntrials,1,1)
    ylabel('spike rate spike/timebin')
    subplot(Ntrials,1,Ntrials)
    xlabel('Time (time bins)')
    hold off
    figure(4)
    cla
    [AX,H1,H2]=plotyy(1:size(Spikes,2),mean(Spikes,1),1:size(Spikes,2),mean(Spikesfilt,1));
    ylabel(AX(1),'Average Spike pattern')
    ylabel(AX(2), 'Filtered Spike pattern')
    xlabel('Time in ms')
    xlim(AX(1),[0 size(Spikes,2)+1])
    xlim(AX(2),[0 size(Spikes,2)+1])
    ylim(AX(1),[0 YMAX2])
    ylim(AX(2),[0 YMAX2])
    set(H1, 'LineWidth', 1.5)
    set(H2, 'LineWidth', 1.5)
    set(AX(1), 'LineWidth', 1.5)
    set(AX(2), 'LineWidth', 1.5)
    title(sprintf('average result of the gaussian filtering over trials with alpha=%d',Alpha_param));
end

end
