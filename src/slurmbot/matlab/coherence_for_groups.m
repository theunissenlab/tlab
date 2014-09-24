function [meanSpikeRate, wholeInfoVals, groupInfoVals] =  coherence_for_groups(cellDir, stimsDir, stimType)

    %% get datafiles in cell directory
    datasets = find_datasets(cellDir, stimsDir, stimType);    
    stimFiles = datasets{1}.srPairs.stimFiles;
    respFiles = datasets{1}.srPairs.respFiles;
    dataDir = datasets{1}.dirname;

    
    %% preprocess the sound to produce spectrograms
    preprocType = 'stft';
    preprocDesc = 'stft.default';
    preprocDir = fullfile(dataDir, 'preproc');
    [s,mess,messid] = mkdir(preprocDir);  

    preprocStimParams = struct;      
    preprocStimParams.tfType = preprocType;
    tfParams = struct;               
    tfParams.high_freq = 8000;       
    tfParams.low_freq = 250;         
    tfParams.log = 1;                
    tfParams.dbnoise = 80;           
    tfParams.refpow = 0;             
    preprocStimParams.tfParams = tfParams;
    preprocStimParams.outputDir = preprocDir;
    preprocStimParams.outputPattern = ['stim.preproc-' preprocDesc '.%d.mat'];
    
    [wholeStim, groupIndex, stimInfo, preprocStimParams] = preprocSound(stimFiles, preprocStimParams);
    
    %% preprocess the responses to produce PSTHs    
    allSpikeTrials = cell(length(respFiles), 1);
    for k = 1:length(respFiles)
        allSpikeTrials{k} = read_spikes_from_file(respFiles{k});
    end

    preprocRespParams = struct;         %create preprocessing struct
    preprocRespParams.units = 'ms';     %our files have units of ms
    preprocRespParams.binSize = 0.001;  %1ms bin size, same sample rate as spectrograms
    preprocRespParams.stimLengths = stimInfo.stimLengths; %needed to compute correct PSTH lengths

    [wholeResponse, groupIndex, respInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);

    
    %% split original spike trials into two PSTHs for purposes of validation
    preprocRespParams.split = 1;
    [wholeSplitResponse, respSplitInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);

    
    numGroups = length(unique(groupIndex));
    avgNumTrials = mean(respInfo.numTrials);
    
    freqCutoff = 90;
    windowSizes = [0.200];
    
    wholeInfoVals = zeros(1, length(windowSizes));
    groupInfoVals = zeros(numGroups, length(windowSizes));
    
    for m = 1:length(windowSizes)
        
        infoWindowSize = windowSizes(m);
    
        cStruct = compute_coherence_bound(wholeSplitResponse(1, :), wholeSplitResponse(2, :), avgNumTrials, stimInfo.sampleRate, freqCutoff, infoWindowSize);        
        wholeInfoVals(m) = cStruct.info;
        clear cStruct;
        
        for k = 1:numGroups
            gIndx = find(groupIndex == k);
            psthHalf1 = wholeSplitResponse(1, gIndx);
            psthHalf2 = wholeSplitResponse(2, gIndx);

            cgStruct = compute_coherence_bound(psthHalf1, psthHalf2, avgNumTrials, stimInfo.sampleRate, freqCutoff, infoWindowSize);        
            groupInfoVals(k, m) = cgStruct.info;
            clear cgStruct;
        end
    
    end
    
    meanSpikeRate = mean(wholeResponse) * stimInfo.sampleRate;
    meanInfoVal = mean(wholeInfoVals);
    meanGroupInfoVals = mean(groupInfoVals, 2);
    
    outputDir = fullfile(dataDir, 'output');
    
    
    [s,mess,messid] = mkdir(outputDir);  
    outputFileName = sprintf('info_vals.txt');
    outputFilePath = fullfile(outputDir, outputFileName);

    fid = fopen(outputFilePath, 'w');
    fprintf(fid, 'meanSpikeRate=%f\n', meanSpikeRate);
    fprintf(fid, 'meanInfoVal=%f\n', meanInfoVal);
    fprintf(fid, 'meanGroupInfoVals=');
    for k = 1:length(meanGroupInfoVals)
        if k ~= 1
            fprintf(fid, ',');
        end
        fprintf(fid, '%f', meanGroupInfoVals(k));
    end
    fprintf(fid, '\n');
    
    fclose(fid);
    
    