function outputs = neural_discrimination(birdname, brainregion, cellname, stimtype)
%[nfiles, stimrate, stdstimrate, backrate, stdbackrate, avgzscore, stdzscore, avgdprime, stddprime, percorrect, avgrank, info, noiseinfo, rateinfo, info1, info2, info2Interval]
% Calculates 3 measures of neural discrimination, the within dprime, the
% percent correct of an ideal oberver based on VR method and the
% information based on Gamma model of spiking neuron
% read all the spike arrival times in a folder

% Read input data from data base
[nfiles spike_times stim_len] = read_all_spikes(birdname, brainregion, cellname, stimtype);

% Window sizes for calculations that depend on window length
winSize = [1 3 5 10 30 50 100];    % Window size for static Info calculation and for ideal observer
% winSize = [50];    % Window size for static Info calculation and for ideal observer
ns = length(winSize);

% Initialize all return values to zero
%***DPrime Calculations***
stimrate = 0;
stdstimrate = 0;
backrate = 0;
stdbackrate = 0;
avgzscore = 0;
stdzscore = 0;
avgdprime = 0;
stddprime = 0;

%***Ideal Observer Calculations (by Trial)***
percorrect = zeros(1,ns);
mi_confusion = zeros(1,ns);

zdistT = zeros(1,ns);               % Uses calculated std for self and other in KL Divergence calc
pzdistT = zeros(1,ns);
mizdistT = zeros(1,ns);

pzdistT2 = zeros(1,ns);             % Uses std for other for both self and other in KL Divergence calc
mizdistT2 = zeros(1,ns);

zdistS = zeros(1,ns);               % Averages distance distribution across song
pzdistS = zeros(1,ns);
mizdistS = zeros(1,ns);

%***Ideal Observer Calculations (by Trial Template)***
percorrectB = zeros(1,ns);
mi_confusionB = zeros(1,ns);

zdistTB = zeros(1,ns);              % Uses std for other for both self and other in KL Divergence calc
pzdistTB = zeros(1,ns);
mizdistTB = zeros(1,ns);
mizdistncTB = zeros(1,ns);

zdistSB = zeros(1,ns);              % Recenters "other" distrubtions to zero to minimize additive variance across songs
pzdistSB = zeros(1,ns);             % and collapses distances across all trials of a single stimulus
mizdistSB = zeros(1,ns);
mizdistncSB = zeros(1,ns);


%***Gamma Model Calculations***
gamma_mutual_info = 0;
gamma_noise_entropy = 0;
gamma_spike_rate = 0;
gamma_const = 0;
rate_bandwidth = 0;
rate_gamma = 0;
fano_factor = 0;

%***Rate Information Calculations***
rate_info_biased = zeros(1,ns);
rate_info_bcorr = zeros(1,ns);
rate_info_stderr = zeros(ns,2);


if nfiles
   % Calculate dprime within
   [stimrate stdstimrate backrate stdbackrate avgzscore stdzscore avgdprime stddprime] = dprime_within(nfiles, spike_times, stim_len);

    % Calculate the ideal oberserver VR metric.
    for is=1:ns
        [pc, mi_conf, zdT, pzdT, pzdT2, mi_zdT, mi_zdT2, zdS, pzdS, mi_zdS] = info_distance(nfiles, spike_times, stim_len, @VR_distance, winSize(is));
        percorrect(is) = pc;
        mi_confusion(is) = mi_conf;
        zdistT(is) = zdT;
        pzdistT(is) = pzdT;
        mizdistT(is) = mi_zdT;
        pzdistT2(is) = pzdT2;
        mizdistT2(is) = mi_zdT2;
        zdistS(is) = zdS;
        pzdistS(is) = pzdS;
        mizdistS(is) = mi_zdS;
    end
    
    
    % With VR_distanceB
    for is=1:ns
        
        % Reinsert this if using the confusion matrix in version B
        [pc, mi_conf, zdT, pzdT, mi_zdT, mi_zdT_nc, zdS, pzdS, mi_zdS, mi_zdS_nc] = info_distanceB(nfiles, spike_times, stim_len, @VR_distanceB, winSize(is));
        percorrectB(is) = pc;
        mi_confusionB(is) = mi_conf;
        
        zdistTB(is) = zdT;
        pzdistTB(is) = pzdT;
        mizdistTB(is) = mi_zdT;
        mizdistncTB(is) = mi_zdT_nc;
        zdistSB(is) = zdS;
        pzdistSB(is) = pzdS;
        mizdistSB(is) = mi_zdS;
        mizdistncSB(is) = mi_zdS_nc;
    end
    

    % Gamma model information calculation
    %onset_time = 200;  % the first 200 ms of the response are removed from the data set
    onset_time = 0;
    spiketrain = spike_times_to_train(nfiles, spike_times, stim_len, onset_time);
    [info, noiseentropy, totalentropy, gamma_const, rate_bandwidth, rate_gamma, fano_factor]= gamma_info(spiketrain);

    gamma_mutual_info = info(1)*1000;                    % Info rates in bits/second
    gamma_noise_entropy = noiseentropy(1)*1000;       % Noiseentropy in bits/second
    gamma_spike_rate = mean(mean(spiketrain))*1000;  % Spiking rate is spikes/second

    % Rate Information
    [rate_info_biased, rate_info_bcorr, rate_info_stderr] = findInfoSR(spiketrain,winSize);

end


outputs = struct('nfiles', nfiles, 'winsize', winSize, 'stimrate', stimrate, 'stdstimrate', stdstimrate, 'backrate', backrate, 'stdbackrate', stdbackrate, ...
    'avgzscore', avgzscore, 'stdzscore', stdzscore, 'avgdprime', avgdprime, 'stddprime', stddprime, ...
    'percorrect', percorrect, 'mi_confusion', mi_confusion, 'zdistT', zdistT, 'pzdistT', pzdistT, 'mizdistT', mizdistT, ...
                                                                              'pzdistT2', pzdistT2, 'mizdistT2', mizdistT2, ...
                                                            'zdistS', zdistS, 'pzdistS', pzdistS, 'mizdistS', mizdistS, ...
    'percorrectB', percorrectB, 'mi_confusionB', mi_confusionB, 'zdistTB', zdistTB, 'pzdistTB', pzdistTB, 'mizdistTB', mizdistTB, 'mizdistncTB', mizdistncTB, ...
                                                                'zdistSB', zdistSB, 'pzdistSB', pzdistSB, 'mizdistSB', mizdistSB, 'mizdistncSB', mizdistncSB, ...                                                                
    'gamma_mutual_info', gamma_mutual_info, 'gamma_noise_entropy', gamma_noise_entropy, 'gamma_spike_rate', gamma_spike_rate,...,
    'gamma_constant', gamma_const, 'rate_bandwidth', rate_bandwidth, 'rate_gamma', rate_gamma, 'fano_factor', fano_factor, ...
    'rate_info_biased',rate_info_biased, 'rate_info_bcorr',rate_info_bcorr,'rate_info_stderr',rate_info_stderr);






return






