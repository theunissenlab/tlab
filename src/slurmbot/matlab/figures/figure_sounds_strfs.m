function figure_sounds_strfs()



    %% make sure strflab functions are in matlab path, in case they aren't already
    strflabDir = get_function_dir('strflab_tutorial_9_Auditory_Example');
    if isempty(strflabDir)
        error('Cannot find strflab directory!');
    end
    addpath(genpath(strflabDir))


    %% get list of stim files
    dataDir = fullfile(strflabDir, 'fakedata', 'auditory');
    audioFiles = get_filenames(dataDir, 'stim[0-9]*.wav', 1);


    %% preprocess the sound to produce spectrograms
    preprocStimParams = struct;      %create preprocessing param structure
    preprocStimParams.tfType = 'stft'; %use short-time FT
    tfParams = struct;               %create time-frequency params
    tfParams.high_freq = 8000;       %specify max freq to analyze
    tfParams.low_freq = 250;         %specify min freq to analyze
    tfParams.log = 1;                %take log of spectrogram
    tfParams.dbnoise = 80;           %cutoff in dB for log spectrogram, ignore anything below this
    tfParams.refpow = 0;             %reference power for log spectrogram, set to zero for max of spectrograms across stimuli
    preprocStimParams.tfParams = tfParams;

    % make a temporary directory to store preprocessed sound files (should be
    %  specific to parameters for preprocSound)
    tempPreprocDir = tempname();    
    [s,mess,messid] = mkdir(tempPreprocDir);
    preprocStimParams.outputDir = tempPreprocDir;

    [wholeStim, groupIndex, stimInfo, preprocStimParams] = preprocSound(audioFiles, preprocStimParams);


    %% preprocess spike times to produce PSTHs
    respFiles = get_filenames(dataDir, 'spike[0-9]*', 1);
    allSpikeTrials = cell(length(respFiles), 1);
    for k = 1:length(respFiles)
        allSpikeTrials{k} = read_spikes_from_file(respFiles{k});
    end

    preprocRespParams = struct;         %create preprocessing struct
    preprocRespParams.units = 'ms';     %our files have units of ms
    preprocRespParams.binSize = 0.001;  %1ms bin size, same sample rate as spectrograms
    preprocRespParams.stimLengths = stimInfo.stimLengths; %needed to compute correct PSTH lengths

    [wholeResponse, groupIndex, respInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);


    
    ng = 11;
    tRng = find(groupIndex == ng);
    stim = wholeStim(tRng, :);
    resp = wholeResponse(tRng);

    tInc = 1 / stimInfo.sampleRate;
    t = 0:tInc:(size(stim, 1)-1)*tInc;

    figure; hold on;
    imagesc(t, stimInfo.f, stim'); axis tight;
    axis xy;
    v_axis = axis;
    v_axis(1) = min(t); v_axis(2) = max(t);
    v_axis(3) = min(stimInfo.f); v_axis(4) = max(stimInfo.f);
    axis(v_axis);
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    
    
    [input, sampleRate, depth] = wavread(fullfile(dataDir, 'stim11.wav'));
    figure; hold on;
    plot(input, 'k-');
    axis tight;
    
    
    figure; hold on;
    plot(resp, 'k-');
    axis tight;
    
    
    ng2 = 5;
    tRng = find(groupIndex == ng2);
    stim = wholeStim(tRng, :);
    resp = wholeResponse(tRng);

    tInc = 1 / stimInfo.sampleRate;
    t = 0:tInc:(size(stim, 1)-1)*tInc;

    figure; hold on;
    imagesc(t, stimInfo.f, stim'); axis tight;
    axis xy;
    v_axis = axis;
    v_axis(1) = min(t); v_axis(2) = max(t);
    v_axis(3) = min(stimInfo.f); v_axis(4) = max(stimInfo.f);
    axis(v_axis);
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    
    figure; hold on;
    plot(resp, 'k-');
    axis tight;
    
    
    
    
    %% display stim and responses
    displayStimResp = 0;

    if displayStimResp
        nGroups = length(unique(groupIndex));
        for k = 1:nGroups

            tRng = find(groupIndex == k);

            stim = wholeStim(tRng, :);
            resp = wholeResponse(tRng);

            tInc = 1 / stimInfo.sampleRate;
            t = 0:tInc:(size(stim, 1)-1)*tInc;

            figure; hold on;

            %plot spectrogram
            subplot(2, 1, 1);
            imagesc(t, stimInfo.f, stim'); axis tight;
            axis xy;
            v_axis = axis;
            v_axis(1) = min(t); v_axis(2) = max(t);
            v_axis(3) = min(stimInfo.f); v_axis(4) = max(stimInfo.f);
            axis(v_axis);
            xlabel('Time (s)'); ylabel('Frequency (Hz)');

            %plot PSTH
            subplot(2, 1, 2);
            plot(t, resp, 'k-'); axis([0 max(t) 0 1]);
            xlabel('Time (s)'); ylabel('P[spike]');

            title(sprintf('Pair %d', k));
        end
    end