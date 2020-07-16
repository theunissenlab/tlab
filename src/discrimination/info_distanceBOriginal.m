function [percorrect, mi_confusion, zdT, pzdT, mi_zdT, zdS, pzdS, mi_zdS] = info_distanceB(nfiles, spike_times, stim_len, distfunc, distfuncparam)
% Calculates the percent correct, the mutual information from the confusion
% matrix, the average z-score, average probalilty and average mutual
% information from the KL-divergence anthropic approximation
% distfunc is the distance function (provided by user), distfuncparam is
% the parameters that go with that function.

debug_flg = 0;                      % Turns debugging on and off

confusion_matrix = zeros(nfiles);   % Stores a matrix of confusiton vectors for each song file
zdT = 0;                                % Stores the updating z-score for a trial
pzdT = 0;                               % Stores the percentage correct based on a gaussian model (by trial)
mi_zdT = 0;                             % Stores the mutual information from the KL-diverence (by trial)
ncountT = 0;                            % Stores the updated number of TRIALS analyzed
zdS = 0;                                % Stores the updating z-score for a song
pzdS = 0;                               % Stores the percentage correct based on a gaussian model (by song)
mi_zdS = 0;                             % Stores the mutual information from the KL-diverence (by song)
ncountS = 0;                            % Stores the updated number of SONGS analyzed

% Compares each trial from each song file to a template created from the
% other trials from the same song, as well as to templates generated from
% combined trials from other songs.
for nf1=1:nfiles
 
    nt1 = length(spike_times{nf1});     % Calculate the number of trials in this file
  
    distance_self = zeros(1, nt1);      % Initialize the self distances for each trial in this file
    distance_other = zeros(1, nfiles-1); % Initialize the other distances for ech trial in this file
    all_other_distance = [];            % This stores all other distances in a single vector
    
    
    
    % Obtain the vector of distances for each spike train from this song to 
    % a template of itself (this song) and the templates created by other songs
    for it1=1:nt1

        % Create a matrix of spike trains from all self trials
        for it=1:nt1
            temp_spikes(it) = spike_times{nf1}(it);
        end

        % Remove the current trial from the collection of spike trains
        temp_spikes(:,it1) = [];
        
        % Concatenate into a single spike train
        temp_spike_times = vertcat(temp_spikes{1:nt1-1});

        % Calculate distance from this spike train to the concatenated
        % spike train of all other "self" trials
        distance_self(it1) = distfunc(spike_times{nf1}{it1}, temp_spike_times, stim_len(nf1), stim_len(nf1), distfuncparam, nt1-1);
        
        % Calculate the distance to templates from the other stimulus
        % presentations
        for nf2=1:nfiles

            % Skip current file
            if nf2 == nf1
                continue;
            end
            
            % Calculate the number of trials in the "other" file
            nt2 = length(spike_times{nf2});
            
            % Initialize the array for recentering the other distributions
            % to zero
            zeroed_other = zeros(1, nt2);       
            
            % Create an empty vector for storing the distances between the
            % current "self" trial and the templates generated from the current "other"
            % file
            distance_other{nf2} = zeros(1, nt2);
            
            % Create a matrix of spike trains from all trials for this
            % stim
            for it=1:nt2
                temp_spikes(it) = spike_times{nf2}(it);
            end
            
%             % Generate a template for each "other" trial, in which all
%             % trials except the current "other" trial are combined and
%             % calculate the distance from the current "self" trial and the
%             % current "other" template
%             for it2=1:nt2             
%                 % Remove the curent "other" trial from the collection of spike trains
%                 temp_spikes2 = temp_spikes;
%                 temp_spikes2(:,it2) = [];
%         
%                 % Concatenate into a single spike train
%                 temp_spike_times = vertcat(temp_spikes2{1:nt2-1});
%                 
%                 % Calculate the distance between the current "self" trial
%                 % and the "other" template
%                 distance_other{nf2}(it2) = distfunc(spike_times{nf1}{it1}, temp_spike_times, stim_len(nf1), stim_len(nf2), distfuncparam, nt2-1);
% 
%             end

            % Generate a template (at random) from each "other" stimulus
            % (using a random 9 of the stimulus' trials)
            remove_trial = round(rand(1)*nt2+0.5);
            
            temp_spikes(:, remove_trial) = [];
                           
            temp_spike_times = vertcat(temp_spikes{1:nt2-1});
                
            % Calculate the distance between the current "self" trial
            % and the "other" template
            distance_other{nf2}(it2) = distfunc(spike_times{nf1}{it1}, temp_spike_times, stim_len(nf1), stim_len(nf2), distfuncparam, nt2-1);

      

            
            % Recenter the distribution of distances to zero
            zeroed_other = distance_other{nf2}-(mean(distance_other{nf2}));

            % Add this distance distribution to the master list of "other" distances
            % for this song
            all_other_distance = [all_other_distance zeroed_other];
        end

        % Calculate the confusion vector for the current "self" trial
        % Compares distance to  "self" template to a random selection of nsamp
        % distances to "other" templates
        confusion_vector = zeros(1, nfiles);
        nsamp = 100;   % Number of ramdom samples from other distribution

        for is=1:nsamp
            distance_samp = zeros(1,nfiles);                % Stores the random distance samples for camparison
            for nf2=1:nfiles
                if nf2 == nf1                               % If on the "self" file...
                    distance_samp(nf2)= distance_self(it1); % Store the "self" distance for the current self trial
                else
                    nt2 = length(spike_times{nf2});         % Calculate the number of trials in the current "other" file
                    it2 = round(rand(1)*nt2+0.5);           % Select a random distance from the "other" file
                    if it2 <= 0
                        fprintf(1, 'Warning: index boundary reached it2 = %d/n', it2);
                        it2 = 1;
                    end
                    if it2 >= nt2+1
                        fprintf(1, 'Warning: index boundary reached it2 = %d/n', it2);
                        it2 = nt2;
                    end
                    distance_samp(nf2) = distance_other{nf2}(it2);  % Store the "other" distance for the current other file
                end
            end

            % Calculate the minimum distance in the randomly selected
            % sample
            min_dist = min(distance_samp);
            
            % Determine which song had the minimum distance
            best_song = find(distance_samp == min_dist);
            
            % Update the confusion vector by increasing the best song's
            % cell by one
            confusion_vector(best_song) = confusion_vector(best_song) + 1;
        end
 
        % Add the current trial's confusion vector to this song's vecotr in the confusion matrix
        confusion_matrix(nf1, :) = confusion_matrix(nf1, :) + confusion_vector;
        
        
         % Calculate the Z scores for this trial
        mean_selfT = distance_self(it1);
        mean_otherT = mean([distance_other{:}]);
        std_otherT = std([distance_other{:}]);
        std_selfT = std_otherT;
                
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
        ncountT = ncountT + 1;                                          % Up the count of analyzed trials
        
        % Debugging figures
        if (debug_flg >= 2)
            figure(1);
            [n_other, x_out] = hist([distance_other{:}], 20);
            bar(x_out, n_other./sum(n_other));
            hold on;
            [n_self, x_out] = hist(distance_self, x_out);
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
            fprintf(1, '\tAverage distance to self = %f\n', mean(distance_self));
            fprintf(1, '\tAverage distance to other = %f\n', mean(mean([distance_other{:}])));
            fprintf(1,'\tStd_other = %f Std_self = %f kld = %f\n', std_other, std_self, kl_divergence./log(2));

            pause();
        end

    end
    
    
    
    
    
    
    % Calculate Z scores for this song
    mean_selfS = mean(distance_self);
    std_selfS = std(distance_self);
    mean_otherS = mean(all_other_distance + mean_selfS);
    std_otherS = std(all_other_distance);
    z_scoreS = (mean_selfS - mean_otherS)./std_otherS;

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
    
    % ***** I believe that this calculation inflates the actual info from
    % kl_divergence... it uses re-centered "other" distance distributions
    % (currently re-centered to zero, and then given an artificial mean of
    % mean_selfS)*****
    kl_divergence = log(std_otherS./std_selfS) + ...                    % Calculate the kl_divergence for this song
        (mean_otherS - mean_selfS)^2./(2*std_otherS.^2) + ...
        (std_selfS.^2 - std_otherS.^2)./(2*std_otherS.^2);
    mi_zdS = mi_zdS + kl_divergence./log(2);                            % Add the mutual info from the kl_divergence for htis song to the
                                                                        % cumulative total -- The log(2) is to go from nits to bits.
    
    ncountS = ncountS + 1;                                              % Up the count of analyzed songs

end



% Normalize the confusion matrix
confusion_matrix = confusion_matrix ./sum(sum(confusion_matrix));   

% Get the mutual information from the confusion matrix
mi_confusion = info_matrix(confusion_matrix);

% Calulate the percent correct by summing the percentage of times the
% calculated song matched the actual song (sum the diagonal ie 1.1, 2.2,
% 3.3)
percorrect = sum(diag(confusion_matrix));

% Calculate the average z-score, guassian percent correct, and
% kl_divergence information
zdT = zdT./ncountT;
pzdT = pzdT./ncountT;
mi_zdT = mi_zdT./ncountT;
zdS = zdS./ncountS;
pzdS = pzdS./ncountS;
mi_zdS = mi_zdS./ncountS;

if (debug_flg >=1 )
    figure(3);
    imagesc(confusion_matrix);
    colorbar;
    xlabel('Actual Song');
    ylabel('Model Song');
    title('Confusion Matrix');
    fprintf(1, 'Pcorrect = %f, MIconfusion = %f\n zd = %f pzd = %f MIzd = %f\n', percorrect, mi_confusion, zd, pzd, mi_zd);
end



return