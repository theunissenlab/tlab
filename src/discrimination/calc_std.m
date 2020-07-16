function [stdSelf stdOther] = calc_std(nfiles, spike_times, stim_len, distfunc, distfuncparam)
% Calculates the percent correct, the mutual information from the confusion
% matrix, the average z-score, average probalilty and average mutual
% information from the KL-divergence anthropomorphic approximation
% distfunc is the distance function (provided by user), distfuncparam is
% the parameters that go with that function.

stdSelf = zeros(nfiles, 10);
stdOther = zeros(nfiles, 10);


% Calculates the distance between each trial from each song file to other 
% trials in the same file and to each trial in all the other files 
for nf1=1:nfiles
    nt1 = length(spike_times{nf1});         % Calculate the number of trials in this file
        

    % Obtain the vector of distances for each spike train with self (this song) and a cell
    % array (to deal with variable number of trials) of distances to other
    % songs
    for it1=1:nt1

        distance_self = zeros(1, nt1-1);        %Stores the distance from this trial to all other trials of the SAME song
        distance_other = cell(1, nfiles-1);     %Stores the distance from this trial to all other trials of all OTHER song


        % Calculate the distance between the current trial and all other
        % trials of the same song
        it_d = 1;
        for it2=1:nt1
            
            %Skip current trial
            if it1 == it2
                continue;
            end
            
            % Add the calculated distance to this trials "self" distances
            distance_self(it_d) = distfunc(spike_times{nf1}{it1}, spike_times{nf1}{it2}, stim_len(nf1), stim_len(nf1), distfuncparam);
            
          
            % Up the count of distances
            it_d = it_d + 1;
        end

        % Calculate the distance from this trial to all trials of all other
        % files
        for nf2=1:nfiles
            
            % Skip current file
            if nf2 == nf1
                continue;
            end
            
            % Calculate the number of trials in the comparision file
            nt2 = length(spike_times{nf2});
            
            % Create an empty vector for storing the distances between the
            % current "self" trial and all trials in the "other" file
            distance_other{nf2} = zeros(1, nt2);
            for it2=1:nt2
                
                % Add the calculated distance to this trial's "other"
                % distances
                distance_other{nf2}(it2) = distfunc(spike_times{nf1}{it1}, spike_times{nf2}{it2}, stim_len(nf1), stim_len(nf2), distfuncparam);
            end
        end

 
     
        
        % Calculate the Z scores for this trial
        stdSelf(nf1, it1) = std(distance_self);
        stdOther(nf1, it1) = std([distance_other{:}]);        
 
      

    end
   
    
end



return