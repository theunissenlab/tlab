function outputs = neural_discrimination_limited(birdname, brainregion, cellname, stimtype, nstims, nperms)
%[nfiles, stimrate, stdstimrate, backrate, stdbackrate, avgzscore, stdzscore, avgdprime, stddprime, percorrect, avgrank, info, noiseinfo, rateinfo, info1, info2, info2Interval]
% Calculates 3 measures of neural discrimination, the within dprime, the
% percent correct of an ideal oberver based on VR method and the
% information based on Gamma model of spiking neuron
% read all the spike arrival times in a folder

% Read input data from data base
[nfiles spike_times stim_len] = read_all_spikes(birdname, brainregion, cellname, stimtype);

if (nstims > nfiles)
    fprintf(1, 'Error: number of possible stims %d is smaller than requested limited %d', nfiles, nstims);
    return;
end


% Window sizes for calculations that depend on window length
winSize = [20];    % Window size for static Info calculation and for ideal observer
% winSize = [50];    % Window size for static Info calculation and for ideal observer
ns = length(winSize);

% Initialize all return values to zero


%***Ideal Observer Calculations (by Trial Template)***
percorrectB = zeros(nperms,ns);
mi_confusionB = zeros(nperms,ns);

zdistTB = zeros(nperms,ns);              % Uses std for other for both self and other in KL Divergence calc
pzdistTB = zeros(nperms,ns);
mizdistTB = zeros(nperms,ns);
mizdistncTB = zeros(nperms,ns);

zdistSB = zeros(nperms,ns);              % Recenters "other" distrubtions to zero to minimize additive variance across songs
pzdistSB = zeros(nperms,ns);             % and collapses distances across all trials of a single stimulus
mizdistSB = zeros(nperms,ns);
mizdistncSB = zeros(nperms,ns);

for iperm=1:nperms
    file_ind = randperm(nfiles);
    spike_times_limited = cell(1,nstims);
    stim_len_limited = zeros(1,nstims);
    for ii=1:nstims
       spike_times_limited(ii) = spike_times(file_ind(ii));
       stim_len_limited(ii) = stim_len(file_ind(ii));
    end
    % With VR_distanceB
    for is=1:ns
        
        % Reinsert this if using the confusion matrix in version B
        [pc, mi_conf, zdT, pzdT, mi_zdT, mi_zdT_nc, zdS, pzdS, mi_zdS, mi_zdS_nc] = info_distanceB(nstims, spike_times_limited, stim_len_limited, @VR_distanceB, winSize(is));
        percorrectB(iperm, is) = pc;
        mi_confusionB(iperm, is) = mi_conf;
        
        zdistTB(iperm, is) = zdT;
        pzdistTB(iperm, is) = pzdT;
        mizdistTB(iperm, is) = mi_zdT;
        mizdistncTB(iperm, is) = mi_zdT_nc;
        zdistSB(iperm, is) = zdS;
        pzdistSB(iperm, is) = pzdS;
        mizdistSB(iperm, is) = mi_zdS;
        mizdistncSB(iperm, is) = mi_zdS_nc;
    end
    

end


outputs = struct('nfiles', nfiles, 'nstims', nstims, 'nperms', nperms, 'winsize', winSize,  ...
      'percorrectB', percorrectB, 'mi_confusionB', mi_confusionB, 'zdistTB', zdistTB, 'pzdistTB', pzdistTB, 'mizdistTB', mizdistTB, 'mizdistncTB', mizdistncTB, ...
      'zdistSB', zdistSB, 'pzdistSB', pzdistSB, 'mizdistSB', mizdistSB, 'mizdistncSB', mizdistncSB );


return