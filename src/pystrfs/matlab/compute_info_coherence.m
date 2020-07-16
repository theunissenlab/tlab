function compute_info_coherence(unitFile, outputFile)

    unit = read_unit_file(unitFile);
        
    for k = 1:length(unit.classes)
    
        className = unit.classes{k};
        resps = unit.class_responses(className);
        
        %% combine all responses within class into cell array
        allSpikeTrials = cell(length(resps), 1);
        stimLengths = zeros(length(resps), 1);
        for m = 1:length(resps)
            allSpikeTrials{k} = resps{m}.spikeTrials;
            stimLengths(m) = resps{m}.trialDuration;
        end

        preprocRespParams = struct;         %create preprocessing struct
        preprocRespParams.units = 's';      %files have units of seconds
        preprocRespParams.binSize = 0.001;  %1ms bin size, same sample rate as spectrograms
        preprocRespParams.stimLengths = stimLengths; %needed to compute correct PSTH lengths
        
        [wholeResponse, groupIndex, respInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);


        %% split original spike trials into two PSTHs for purposes of validation
        preprocRespParams.split = 1;
        [wholeSplitResponse, respSplitInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);

        avgNumTrials = mean(respInfo.numTrials);

        windowSizes = [0.200];

        wholeInfoVals = zeros(1, length(windowSizes));

        for m = 1:length(windowSizes)

            infoWindowSize = windowSizes(m);

            cStruct = compute_coherence_bound(wholeSplitResponse(1, :), wholeSplitResponse(2, :), avgNumTrials, stimInfo.sampleRate, freqCutoff, infoWindowSize);        
            wholeInfoVals(m) = cStruct.info;
            clear cStruct;

        end

        meanSpikeRate = mean(wholeResponse) * stimInfo.sampleRate;
        meanInfoVal = mean(wholeInfoVals);
        
    
        
        
        

    end