function outputs = compute_discrimination_unit(h5Path)
    
    %[nstim, stimrate, stdstimrate, backrate, stdbackrate, avgzscore, stdzscore, avgdprime, stddprime, percorrect, avgrank, info, noiseinfo, rateinfo, info1, info2, info2Interval]
    % Calculates 3 measures of neural discrimination, the within dprime, the
    % percent correct of an ideal oberver based on VR method and the
    % information based on Gamma model of spiking neuron
    % read all the spike arrival times in a folder

    %% Read input data from data base
    unit = read_unit_file(h5Path);

    %% Select the protocol (SelX, Callx, Maskx, STRFx...)
    %number of different protocols run for that unit, identify if the one
    %you are asking for is existing in this unit and selecting this one in
    %responses
    nclasses = length(unit.classes);
    
    %% write stuff out to a file
    [path, fname, ext] = fileparts(h5Path);
    outputFile = fullfile(path, 'discrimination_info.h5');
    
    h5 = h5utils();
    f = h5.create(outputFile);    
    for ic=1:nclasses   
        stimClass = unit.classes{ic};
        compute_discrimination_unit_oneclass(unit, stimClass, f);   
    end
    
    h5.close(f);   

end

function compute_discrimination_unit_oneclass(unit, stimType, f)

    responses = unit.class_responses.(stimType);

    %% Construct the input we need for all the following code
    % This is the number of sound files played
    nstim = length(responses);

    %This is the vector of the length of the different stims
    stim_len = zeros (nstim,1);
    for nf=1:nstim
        stim_len(nf) = 1000*responses{nf}.trialDuration;%here stim_len is in ms
    end

    % This is the cell array of the time of spike arrival
    spike_times = cell(1, nstim);
    for nf = 1:nstim
        Trials = responses{nf}.spikeTrials;
        nt = length(Trials);
        spike_times{nf} = cell(1, nt);
        for it = 1:nt        
            spike_times{nf}{it} = Trials{it}*1000; %has to be in ms
        end
    end

    %% Set up parameters for all the calculations    
    % Window sizes for calculations that depend on window length
    winSize = [1 3 5 10 30 50 100];    % Window size for static Info calculation and for ideal observer
    % winSize = [50];    % Window size for static Info calculation and for ideal observer
    ns = length(winSize);

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


    %% Calculate dprime within
    [stimrate stdstimrate backrate stdbackrate avgzscore stdzscore avgdprime stddprime] = dprime_within(nstim, spike_times, stim_len);

    %% Calculate the ideal oberserver VR metric with VR_distanceB
    for is=1:ns

        % Reinsert this if using the confusion matrix in version B
        [pc, mi_conf, zdT, pzdT, mi_zdT, mi_zdT_nc, zdS, pzdS, mi_zdS, mi_zdS_nc, confusion_matrix] = info_distanceB(nstim, spike_times, stim_len, @VR_distanceB, winSize(is));
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


    %% Compute Gamma model information
    %onset_time = 200;  % the first 200 ms of the response are removed from the data set
    onset_time = 0;
    spiketrain = spike_times_to_train(nstim, spike_times, stim_len, onset_time);
    [info, noiseentropy, totalentropy, gamma_const, rate_bandwidth, rate_gamma, fano_factor]= gamma_info(spiketrain);

    gamma_mutual_info = info(1)*1000;                    % Info rates in bits/second
    gamma_noise_entropy = noiseentropy(1)*1000;       % Noiseentropy in bits/second
    gamma_spike_rate = mean(mean(spiketrain))*1000;  % Spiking rate is spikes/second

    %% Compute Rate Information
    [rate_info_biased, rate_info_bcorr, rate_info_stderr] = findInfoSR(spiketrain,winSize);

    %% write data to an hdf5 file (already opened)
    h5 = h5utils();
    stimGrp = sprintf('/%s', stimType);
    h5.create_group(f, stimGrp);

    h5.set_ds(f, stimGrp, 'win_size', winSize);

    rateGrp = sprintf('%s/rate_info', stimGrp);
    h5.create_group(f, rateGrp);
    h5.set_ds(f, rateGrp, 'info_biased', rate_info_biased);
    h5.set_ds(f, rateGrp, 'info_bcorr', rate_info_bcorr);
    h5.set_ds(f, rateGrp, 'info_stderr', rate_info_stderr);

    gammaGrp = sprintf('%s/gamma_info', stimGrp);
    h5.create_group(f, gammaGrp);
    h5.set_ds(f, gammaGrp, 'gamma_mutual_info', gamma_mutual_info);
    h5.set_ds(f, gammaGrp, 'gamma_noise_entropy', gamma_noise_entropy);
    h5.set_ds(f, gammaGrp, 'fano_factor', fano_factor);
    h5.set_ds(f, gammaGrp, 'rate_gamma', rate_gamma);
    h5.set_ds(f, gammaGrp, 'noise_entropy', noiseentropy);
    h5.set_ds(f, gammaGrp, 'total_entropy', totalentropy);
    h5.set_ds(f, gammaGrp, 'info', info);
    h5.set_ds(f, gammaGrp, 'rate_bandwidth', rate_bandwidth);

    catGrp = sprintf('%s/categorical_info', stimGrp);
    h5.create_group(f, catGrp);
    h5.set_ds(f, catGrp, 'percent_correct', percorrectB);
    h5.set_ds(f, catGrp, 'mi_confusion', mi_confusionB);
    h5.set_ds(f, catGrp, 'zdistT', zdistTB);
    h5.set_ds(f, catGrp, 'p_zdistT', pzdistTB);
    h5.set_ds(f, catGrp, 'mi_zdistT', mizdistTB);
    h5.set_ds(f, catGrp, 'zdistS', zdistSB);
    h5.set_ds(f, catGrp, 'p_zdistS', pzdistSB);
    h5.set_ds(f, catGrp, 'mi_zdistS', mizdistSB);
    h5.set_ds(f, catGrp, 'confusion_matrix', confusion_matrix);
end
