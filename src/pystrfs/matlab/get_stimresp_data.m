function stimRespData = get_stimresp_data(unitFile, preprocFile, stimClass)

    h5 = h5utils();
        
    %% get unit data, including responses and stimuli
    ufid = h5.open(unitFile);
    pfid = h5.open(preprocFile);
    
    stimGroup = sprintf('/class_responses/%s', stimClass);
    numStims = h5.get_attr(ufid, stimGroup, 'num_stims');
    
    md5s = cell(numStims, 1);
    allSpikeTrials = cell(numStims, 1);
    stimLengths = zeros(numStims, 1);
    actualStimLengths = zeros(numStims, 1);
    stims = cell(numStims, 1);
    
    totalLength = 0;
    numChannels = -1;
    lowFreq = -1;
    highFreq = -1;
    sampleRate = -1;
    
    allNumTrials = zeros(length(stims), 1);
    
    
    [rootDir, pname, ext] = fileparts(preprocFile);
    is_neurogram = strcmp('neurogram', pname);
    
    if is_neurogram
        h5.close(pfid);
        sampleRate = 1000;
    end
    
    for k = 1:length(stims)
       
        %% get metadata
        respGroup = sprintf('%s/%d', stimGroup, k);
        md5s{k} = h5.get_attr(ufid, respGroup, 'stim_md5');
        stimLengths(k) = h5.get_attr(ufid, respGroup, 'trial_duration');
        numTrials = h5.get_attr(ufid, respGroup, 'num_trials');
        allNumTrials(k) = numTrials;
        
        %% get spike trials
        spikeTrials = cell(numTrials, 1);
        for m = 1:int8(numTrials)
            trialDs = sprintf('%s/%d', respGroup, m);
            spikeTimes = h5.get_ds(ufid, trialDs);        
            spikeTrials{m} = spikeTimes;
        end        
        allSpikeTrials{k} = spikeTrials;        
        
        %% get spectrogram
        if ~is_neurogram
            specGroup = sprintf('/%s', md5s{k});
            stimDs = sprintf('%s/spectrogram', specGroup);        

            lowFreq = h5.get_attr(pfid, specGroup, 'low_freq');
            highFreq = h5.get_attr(pfid, specGroup, 'high_freq');
            sampleRate = h5.get_attr(pfid, specGroup, 'sample_rate');

            stims{k} = h5.get_ds(pfid, stimDs);

            actualStimLengths(k) = size(stims{k}, 2) / sampleRate;
            totalLength = totalLength + size(stims{k}, 2);
            numChannels = size(stims{k}, 1);        
        end
    end
    
    if ~is_neurogram
        h5.close(pfid);
    end
    
    if is_neurogram
        
        [rootDir, ustr, ext] = fileparts(unitFile);
        [rootDir, unitName, ext] = fileparts(rootDir);
        
        allUnits = get_good_units();
        %remove this unit
        unitsToUse = {};
        for k = 1:length(allUnits)
            u = allUnits{k};
            if ~strcmp(u, unitName)
                unitsToUse{end+1} = allUnits{k};
            end            
        end
        
        [wholeStim, groupIndex, ngramInfo] = preproc_neurogram(preprocFile, md5s, unitsToUse, 'best');        
        sz = size(wholeStim);
        numChannels = sz(1);
        actualStimLengths = ngramInfo.stimLengths / sampleRate;
        
        lowFreq = 1;
        highFreq = length(unitsToUse);
    else
        %% concatenate spectrograms
        wholeStim = zeros(numChannels, totalLength);
        cindx = 1;
        for k = 1:length(stims)        
            slen = size(stims{k}, 2);
            eindx = cindx+slen-1;
            wholeStim(:, cindx:eindx) = stims{k};
            cindx = eindx+1;
        end
    end
    
    h5.close(ufid);    
    
    
    
    %% reconstruct frequencies
    freqInc = (highFreq - lowFreq) / numChannels;
    f = lowFreq:freqInc:highFreq;
    
    %% create PSTHs
    preprocRespParams = struct;         %create preprocessing struct
    preprocRespParams.units = 's';      %files have units of s
    preprocRespParams.binSize = 0.001;  %1ms bin size, same sample rate as spectrograms
    preprocRespParams.stimLengths = actualStimLengths; %needed to compute correct PSTH lengths

    [wholeResponse, groupIndex, respInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);

    preprocRespParams.split = 1;
    [wholeSplitResponse, respSplitInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);

    
    %% create output structure
    stimRespData = struct;
    stimRespData.f = f;
    stimRespData.numChannels = numChannels;
    stimRespData.sampleRate = sampleRate;
    stimRespData.allSpikeTrials = allSpikeTrials;
    stimRespData.wholeStim = wholeStim';
    stimRespData.groupIndex = groupIndex;
    stimRespData.wholeResp = wholeResponse;
    stimRespData.wholeSplitResp = wholeSplitResponse;
    stimRespData.numTrials = allNumTrials;
    
    