function [percorrect, mi_confusion, zdT, pzdT, mi_zdT, mi_zdT_nc, zdS, pzdS, mi_zdS, mi_zdS_nc, confusion_matrix, ncountT,Distance_allOtherStims, Response_allStims] = info_distanceB_slvg(nfiles, spike_times, stim_len, distfunc, distfuncparam , debug_flg)

% Reinsert this if not using the confusion matrix!!!

% function [percorrect, mi_confusion, zdT, pzdT, mi_zdT, zdS, pzdS, mi_zdS] = info_distanceB(nfiles, spike_times, stim_len, distfunc, distfuncparam)


% Calculates the percent correct, the mutual information from the confusion
% matrix, the average z-score, average probalilty and average mutual
% information from the KL-divergence anthropic approximation
% distfunc is the distance function (provided by user), distfuncparam is
% the parameters that go with that function.
if nargin<6
    debug_flg = 0;                      % Turns debugging on and off
end

% Find the minimum stim length
templateSize = min(stim_len);

if nargin<5
   distfuncparam = templateSize; % Fix the Window parameter of SR_distanceB or VR_distanceB to the size of the template if nothing is specified 
   s=functions(distfunc);
   if strcmp(s.function, 'VR_distanceB')
        fprintf('WARNING: No Window size was specified, the length of the stimulus is used as a window size for VR_distanceB/n');
   end
end
% Initialize a matrix that will contain the average distance (over trials) of each stim
% (rows) to the other stims (columns)
Distance_allOtherStims = nan(nfiles);
Response_allStims = cell(nfiles,1);

confusion_matrix = zeros(nfiles);   % Stores a matrix of confusiton vectors for each song file
zdT = 0;                                % Stores the updating z-score for a trial
pzdT = 0;                               % Stores the percentage correct based on a gaussian model (by trial)
mi_zdT = 0;                             % Stores the mutual information from the KL-diverence (by trial)
mi_zdT_nc = 0;                          % Stores the mutual information from the KL-divergence without anthropic correction (by trial)
ncountT = 0;                            % Stores the updated number of TRIALS analyzed
zdS = 0;                                % Stores the updating z-score for a song
pzdS = 0;                               % Stores the percentage correct based on a gaussian model (by song)
mi_zdS = 0;                             % Stores the mutual information from the KL-diverence (by song)
mi_zdS_nc = 0;                          % Stores the mutual information from the KL-diverence (by song) without anthropic correction
ncountS = 0;                            % Stores the updated number of SONGS analyzed




% Set the parameters of the distance calculations
s=functions(distfunc);
if distfuncparam==templateSize || unique(round(stim_len))==round(templateSize) %no sliding if the size of the window is of the size of the spike train or for equally long stims
        shiftflg=0;
else
        shiftflg=1;
end
gaussflg=0; % We go for exponential smoothing as a default.
normflg=0; % We go for no normalization of spike patterns by default, so we take average spike rate differences into account

if strcmp(s.function, 'VR_distanceB')
    fprintf(1, 'Info_distanceB is calculating the Van Rossum distance between spike patterns\n');
    if (gaussflg == 1)
       fprintf(1, 'Info_distanceB is applying a gaussian smoothing of winsize %d on spike trains\n', distfuncparam);
    else
       fprintf(1,'Info_distanceB is applying a exponential smoothing of winsize %d on spike trains\n', distfuncparam);
    end
    
    if (normflg == 1)
       fprintf(1, 'Info_distanceB is using distances between templates normalized by rate\n');
    else
       fprintf(1,'Info_distanceB is using distances between non-rate-normalized templates\n');
    end 

    if (shiftflg == 1)
       fprintf(1, 'Info_distanceB is calculating best distances between templates by sliding templates\n');
    else
       fprintf(1,'Info_distanceB is calculating the distance between aligned templates\n');
    end
elseif strcmp(s.function, 'SR_distanceB')
    fprintf(1, 'Info_distanceB is calculating the distance between the spike counts using SR_distanceB\n');
end

%% Get weights for each neuron
nweights = cell_weights(spike_times,stim_len,distfuncparam);%%What is this code doing in case of h5 file only 1 unit per file, so nweights=1
%nweights = ones(1,size(nweights));
% nweights = 1;

% Compares each trial from each song file to a template created from the
% other trials from the same song, as well as to templates generated from
% combined trials from other songs.
for nf1=1:nfiles
    %fprintf(1,'WinSize=%d calculus on file %d out of %d\n', distfuncparam, nf1, nfiles);
 
    nt1 = length(spike_times{nf1});     % Calculate the number of trials in this file
  
    distance_self = zeros(1, nt1);      % Initialize the self distances for each trial in this file
    distance_other = zeros(nt1, nfiles);   % Initialize the other distances for all trials in this file
    response_self = cell(nt1,1);              % Initialize the cell containing the wavforms (VR_distanceB) or spike count (SR_distanceB) for each trial in this file
    all_self_distance = [];                 %Stores all SELF distances for all trials of THIS song
    all_other_distance = [];                %Stores all OTHER distances for all trials of THIS song

    
    
    % Obtain the vector of distances for each spike train from this song to 
    % a template of itself (this song) and the templates created by other songs
    for it1=1:nt1
        %fprintf(1,'Winsize=%d calculus of self and other distances on trial %d out of %d\n', distfuncparam, it1, nt1);

        % Create a matrix of spike trains from all self trials
        temp_spikes=cell(1,nt1);
        for it=1:nt1
            temp_spikes(it) = spike_times{nf1}(it);
        end

        % Remove the current trial from the collection of spike trains
        temp_spikes(:,it1) = [];
        
        % Concatenate into a single spike train
        SizeSpike=1;
        ss=0;
        while SizeSpike<2 && ss<(nt1-1)
            ss=ss+1;
            SizeSpike=length(temp_spikes{ss});
        end
            
        if size(temp_spikes{ss},1)==1 %row vector
            temp_spike_times = horzcat(temp_spikes{1:nt1-1});
            temp_spike_times = temp_spike_times';%convert to column vector
        else  %column vector or only one element or none in each cell
            temp_spike_times = vertcat(temp_spikes{1:nt1-1});
        end
        
        % check dimension of the tested spike train
        Ispike_times=spike_times{nf1}{it1};
        if size(Ispike_times,1)==1
            Ispike_times=Ispike_times';
        end

        % Calculate distance from this spike train to the concatenated
        % spike train of all other "self" trials
        [distance_self(it1),response_self{it1}] = distfunc(Ispike_times, temp_spike_times, stim_len(nf1), stim_len(nf1), distfuncparam, nt1-1, nweights, templateSize, normflg, shiftflg, gaussflg);
        
        
        % Calculate the distance from "self" to a random template from the other stimulus
        % presentations
        for nf2=1:nfiles
            %fprintf(1,'calculating distance other trial %d with file %d\n', it1, nf2);

            % Skip current file
            if nf2 == nf1
                continue;
            end
            
            % Calculate the number of trials in the "other" file
            nt2 = length(spike_times{nf2});
            
            % Initialize the array for recentering the distributions
            % to zero
            zeroed_other = zeros(1, nt2);       
            zeroed_self = 0;
             
            % Create a matrix of spike trains from all trials for this
            % stim
            temp_spikes=cell(1,nt2);
            for it=1:nt2
                temp_spikes(it) = spike_times{nf2}(it);
            end
            
            % Generate a template for a random "other" trial, selected "other" trial are combined and
            % calculate the distance from the current "self" trial and the
            % current "other" template
            it2 = round(rand(1)*nt2+0.5);
            
            % Remove the curent "other" trial from the collection of spike trains
            temp_spikes2 = temp_spikes;
            temp_spikes2(:,it2) = [];
            
            %determine the shape of the vector for the first trial where
            %there is more than a single spike and creectly shape the
            %vector for VR_distanceB
            SizeSpike=1;
            ss=0;
            while SizeSpike<2 && ss<(nt2-1)
                ss=ss+1;
                SizeSpike=length(temp_spikes2{ss});
            end
        
            if size(temp_spikes2{ss},1)==1 %row vector
                temp_spike_times = horzcat(temp_spikes2{1:nt2-1});
                temp_spike_times = temp_spike_times';%convert to column vector
            else %column vector or only one element or none in each cell
                temp_spike_times = vertcat(temp_spikes2{1:nt2-1});
            end
                
            % Calculate the distance between the current "self" trial
            % and the "other" template
            distance_other(nt1,nf2) = distfunc(Ispike_times, temp_spike_times, stim_len(nf1), stim_len(nf2), distfuncparam, nt2-1, nweights,templateSize, normflg, shiftflg, gaussflg);
        end
        
        % distance_other also includes a zero distance for nf1==nf2
        distance_other_trunc = distance_other(nt1,:);
        distance_other_trunc(nf1) = [];
        
        % Recenter the distribution of distances to zero
        mean_distance_other = mean(distance_other_trunc);
        
        zeroed_other = distance_other_trunc-mean_distance_other;
        zeroed_self = distance_self(it1)-mean_distance_other;


        % Add this distance distribution to the master list of "other" distances
        % for this song
        all_other_distance = [all_other_distance zeroed_other];
        all_self_distance = [all_self_distance zeroed_self];
        
        
        
        % Calculate the confusion vector for the current "self" trial
        % Compares distance to  "self" template to a random selection of nsamp
        % distances to "other" templates
        confusion_vector = zeros(1, nfiles);
       
        distance_samp = zeros(1,nfiles);                % Stores the random distance samples for camparison
        for nf2=1:nfiles
            if nf2 == nf1                               % If on the "self" file...
                distance_samp(nf2)= distance_self(it1); % Store the "self" distance for the current self trial
            else
                distance_samp(nf2) = distance_other(nt1,nf2);  % Store the "other" distance for the current other file
            end
        end

        % Calculate the minimum distance in the randomly selected
        % sample
        min_dist = min(distance_samp);

        % Determine which song had the minimum distance
        best_song = find(distance_samp == min_dist);
        bb = length(best_song);

        % Update the confusion vector by increasing the best song's
        % cell by one
        confusion_vector(best_song) = confusion_vector(best_song) + 1/bb;
       
 
        % Add the current trial's confusion vector to this song's vecotr in the confusion matrix
        confusion_matrix(nf1, :) = confusion_matrix(nf1, :) + confusion_vector;
        
        
         % Calculate the Z scores for this trial
        mean_selfT = distance_self(it1);
        mean_otherT = mean_distance_other;
        std_otherT = std(distance_other_trunc);
        std_selfT = std_otherT;
        if (std_otherT == 0)
            fprintf(1, 'Stim %d Trial %d: \n', nf1, it1);
            fprintf(1, 'Critical Warning in info_distance B: trial is equally distant to templates to all other stimuli');
        end
        mean_allT = mean(distance_samp);
        std_allT = std(distance_samp);
                
        z_scoreT = (mean_selfT - mean_otherT)./std_otherT;
      
        % Calculate the percent correct based on the gaussian model
        nsamp = 1000;
        p_valueT = 0;
        
        % Compare nsamp # of randomly generated "self distances" to nsamp # of
        % randomly gerated "other distances" and calculate the percentage of
        % times that the "self distance" is the shortest.        
        for is=1:nsamp
            dsamp_self = randn(1).*std_selfT + mean_selfT;              % Generate a point in the self distribution
            dsamp_other = randn(1,nfiles-1).*std_otherT + mean_otherT;  % Generate a point in the other distribution for each other song
            if (dsamp_self < min(dsamp_other))                          % Determine if the self distance is the shortest
                p_valueT = p_valueT +1;
            end
        end
        p_valueT = p_valueT./nsamp;                                     % Calculate the number of times the self distnace was the shortest
        zdT = zdT + z_scoreT;                                           % Add the current trial z-score to the cumulative total
        pzdT = pzdT + p_valueT;                                         % Add the current (trial) gaussian percentage to the cumulative total
        kl_divergence = log(std_otherT./std_selfT) + ...                % Calculate the kl_divergence for this trial
        (mean_otherT - mean_selfT)^2./(2*std_otherT.^2) + ...
        (std_selfT.^2 - std_otherT.^2)./(2*std_otherT.^2);
        mi_zdT = mi_zdT + kl_divergence./log(2);                        % Add the mutual info from the kl_divergence for htis song to the
                                                                        % cumulative total -- The log(2) is to go from nits to bits.
        % Repeat without anthropic correction
        kl_all =  (mean_allT - mean_selfT)^2./(2*std_allT.^2);
        mi_zdT_nc = mi_zdT_nc + kl_all./log(2);                        % Add the mutual info from the kl_divergence for htis song to the
                                                                        % cumulative total -- The log(2) is to go from nits to bits.
        ncountT = ncountT + 1;                                          % Up the count of analyzed trials
        
        % Debugging figures
        if (debug_flg >= 3)
            figure(1);
            [n_other, x_out] = hist(distance_other_trunc, 20);
            bar(x_out, n_other./sum(n_other));
            hold on;
            [n_self, x_out] = hist(distance_self(it1), x_out);
            bar(x_out, n_self./sum(n_self),'r' );
            legend('Other', 'Own');
            xlabel('VR Distance');
            ylabel('Frequency');
            hold off;

            figure(2);
            bar(confusion_vector./sum(confusion_vector));
            xlabel('Song Number');
            ylabel('Probability');
            title('Confusion vector');


            fprintf(1, 'Song %d Trial %d: \n', nf1, it1);
            fprintf(1, '\tAverage distance to self = %f\n', mean_selfT);
            fprintf(1, '\tAverage distance to other = %f\n', mean_otherT);
            fprintf(1,'\tStd_other = %f Std_self = %f kld = %f kld(nc)=%f\n', std_otherT, std_selfT, kl_divergence./log(2), kl_all./log(2));


            pause();
        end

    end
    
    % Save the average distances (over trials) of that stim to other stims
    Distance_allOtherStims(nf1,:)= mean(distance_other,1);
    
    % Save the average spike pattern or spike count for that stim
    Response_allStims{nf1} = mean(cell2mat(response_self),1);
    
    % Calculate Z scores for this song
    mean_selfS = mean(all_self_distance);
    std_selfS = std(all_self_distance);
    mean_otherS = mean(all_other_distance);
    std_otherS = std(all_other_distance);
    z_scoreS = (mean_selfS - mean_otherS)./std_otherS;
    if (std_selfS == 0)
        if (std_otherS == 0)
            fprintf(1, 'Warning in info_distanceB: all distances to Self for stim %d are equal setting std_selfS to std_otherS\n', nf1);
            std_selfS = std_otherS;
        end
    end
    if (std_otherS == 0)
        fprintf(1, 'Critical Warning in info_distance B: there is no differential responses to stimuli for stim %d', nf1);
    end

    % Calculate the percent correct based on the gaussian model
    nsamp = 1000;
    p_valueS = 0;
    
    % Compare nsamp # of randomly generated "self distances" to nsamp # of
    % randomly gerated "other distances" and calculate the percentage of
    % times that the "self distance" is the shortest.
    for is=1:nsamp
        dsamp_self = randn(1).*std_selfS + mean_selfS;                  % Generate a point in the self distribution
        dsamp_other = randn(1,nfiles-1).*std_otherS + mean_otherS;      % Generate a point in the other distribution for each other song
        if (dsamp_self < min(dsamp_other))                              % Determine if the self distance is the shortest
            p_valueS = p_valueS +1;
        end
    end
    p_valueS = p_valueS./nsamp;                                         % Calculate the number of times the self distnace was the shortest
   
    zdS = zdS + z_scoreS;                                               % Add the current z-score to the cumulative total
    pzdS = pzdS + p_valueS;                                             % Add the current gaussian percentage to the cumulative total
    
    kl_divergence = log(std_otherS./std_selfS) + ...                    % Calculate the kl_divergence for this song
        (mean_otherS - mean_selfS)^2./(2*std_otherS.^2) + ...
        (std_selfS.^2 - std_otherS.^2)./(2*std_otherS.^2);
    mi_zdS = mi_zdS + kl_divergence./log(2);                            % Add the mutual info from the kl_divergence for htis song to the
                                                                        % cumulative total -- The log(2) is to go from nits to bits.
    
    % Repeat kl_divergence without anthropic correction
    all_distances = [all_self_distance all_other_distance];
    mean_allS = mean(all_distances);
    std_allS = std(all_distances);
    kl_all = log(std_allS./std_selfS) + ...                    % Calculate the kl_divergence for this song
        (mean_allS - mean_selfS)^2./(2*std_allS.^2) + ...
        (std_selfS.^2 - std_allS.^2)./(2*std_allS.^2);
    mi_zdS_nc = mi_zdS_nc + kl_all./log(2);                            % Add the mutual info from the kl_divergence for htis song to the
                                                                        % cumulative total -- The log(2) is to go from nits to bits.                                                                    
                                                                        
    ncountS = ncountS + 1;                                              % Up the count of analyzed songs
            % Debugging figures
        if (debug_flg >= 2)
            figure(3);
            x_out = min(all_distances):0.01:max(all_distances)
            [n_other, x_out] = hist(all_other_distance, x_out);
            bar(x_out, n_other./sum(n_other));
            hold on;
            dx = x_out(2)-x_out(1);
            [n_self, x_out] = hist(all_self_distance, x_out);
            bar(x_out, n_self./sum(n_self),'r' );
            plot(x_out, normpdf(x_out, mean_selfS, std_selfS).*dx, 'k', 'LineWidth',2);
            plot(x_out, normpdf(x_out, mean_otherS, std_otherS).*dx, 'k', 'LineWidth',2);
            legend('Other', 'Self');
            xlabel('Zerod VR Distance');
            ylabel('Frequency');
            title(sprintf('KL divergence %.2f bits\n', kl_divergence./log(2)));
            hold off;
            fprintf(1, 'Song %d\n', nf1);
            fprintf(1, '\tAverage distance to self = %f\n', mean_selfS);
            fprintf(1, '\tAverage distance to other = %f\n', mean_otherS);
            fprintf(1,'\tStd_other = %f Std_self = %f kld = %f kld(no correction) = %f\n', std_otherS, std_selfS, kl_divergence./log(2), kl_all./log(2));
            pause();
        end

end

if ncountT~=round(sum(sum(confusion_matrix)))
    fprintf('WARNING MAJOR ERROR of info_distanceB in the number of trials. There are %d trials in spike_times and %d trials used to construct the confusion matrix\n', ncountT, sum(sum(confusion_matrix)));
end

% Normalize the confusion matrix
confusion_matrix = confusion_matrix ./sum(sum(confusion_matrix));

% Get the mutual information from the confusion matrix
mi_confusion = info_matrix(confusion_matrix);

% Calulate the percent correct by summing the percentage of times the
% calculated song matched the actual song (sum the diagonal ie 1.1, 2.2,
% 3.3)
percorrect = sum(diag(confusion_matrix));

% Calculate the average z-score, gaussian percent correct, and
% kl_divergence information
zdT = zdT./ncountT;
pzdT = pzdT./ncountT;
mi_zdT = mi_zdT./ncountT;
mi_zdT_nc = mi_zdT_nc./ncountT;
zdS = zdS./ncountS;
pzdS = pzdS./ncountS;
mi_zdS = mi_zdS./ncountS;
mi_zdS_nc = mi_zdS_nc./ncountS;


if (debug_flg >=1 )
    figure(4);
    imagesc(confusion_matrix);
    colorbar;
    xlabel('Model Song');
    ylabel('Actual Song');
    title('Confusion Matrix');
    fprintf(1, 'Pcorrect = %f, MIconfusion = %f\n Song Calculation zd = %f pzd = %f MIzd = %f MIzd(nc) = %f\n Trial Calculation zd = %f pzd = %f MIzd = %f MIzd(nc) = %f\n', ...
        percorrect, mi_confusion, zdS, pzdS, mi_zdS, mi_zdS_nc, zdT, pzdT, mi_zdT, mi_zdT_nc);
end



return