function preprocData = preproc_stim_response(cellDir, stimsDir, zscore, stimTypes, tfType, tfParams, preprocDesc)

    if ~isstruct(zscore)
        zs = zscore;
        zscore = struct;
        zscore.stim = zs;
        zscore.resp = zs;
    end

    %% get datafiles in cell directory
    datasets = find_datasets(cellDir, stimsDir, stimTypes);    
    stimFiles = datasets{1}.srPairs.stimFiles;
    respFiles = datasets{1}.srPairs.respFiles;
    dataDir = datasets{1}.dirname;
   
    
    %% preprocess the sound to produce spectrograms
    preprocType = tfType;
    preprocDir = fullfile(dataDir, 'preproc');
    [s,mess,messid] = mkdir(preprocDir);  

    preprocStimParams = struct;      
    preprocStimParams.tfType = preprocType;
    preprocStimParams.tfParams = tfParams;
    preprocStimParams.outputDir = preprocDir;
    preprocStimParams.overwrite = 1;
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

    
    %% z-score the stimulus and response
    stimInfo.zscored = zscore.stim;
    respInfo.zscored = zscore.resp;
    if zscore.stim
        fprintf('z-scoring stimulus...\n');
        wholeStim = norm_std_mean(wholeStim);
    end
    if zscore.resp
        fprintf('z-scoring response...\n');
        respInfo.mean = mean(wholeResponse);
        respInfo.std = std(wholeResponse);
        wholeResponse = (wholeResponse - respInfo.mean) / respInfo.std;
    end
    
    
    %% split original spike trials into two PSTHs for purposes of validation
    preprocRespParams.split = 1;
    [wholeSplitResponse, respSplitInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);

    
    preprocData.stim = wholeStim;
    preprocData.resp = wholeResponse;
    preprocData.splitResp = wholeSplitResponse;
    preprocData.groupIndex = groupIndex;
    preprocData.stimParams = preprocStimParams;
    preprocData.stimInfo = stimInfo;
    preprocData.respParams = preprocRespParams;
    preprocData.respInfo = respInfo;
    preprocData.dataDir = dataDir;
    