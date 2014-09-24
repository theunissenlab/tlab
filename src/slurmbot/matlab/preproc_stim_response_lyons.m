function preprocData = preproc_stim_response_lyons(cellDir, stimsDir, zscore, stimTypes)

    if nargin < 3
        zscore = 1;
    end
    
    if nargin < 4
        stimTypes = 'conspecific|zfsongs';
    end

    %% get datafiles in cell directory
    datasets = find_datasets(cellDir, stimsDir, stimTypes);    
    stimFiles = datasets{1}.srPairs.stimFiles;
    respFiles = datasets{1}.srPairs.respFiles;
    dataDir = datasets{1}.dirname;
   
    
    %% preprocess the sound to produce spectrograms
    preprocType = 'lyons';
    preprocDir = fullfile(dataDir, 'preproc');
    [s,mess,messid] = mkdir(preprocDir);  

    preprocStimParams = struct;    
    preprocStimParams.overwrite = 1;
    preprocStimParams.tfType = preprocType;
    
    tfParams = struct;  
    tfParams.low_freq = 300;
    tfParams.high_freq = 15000;
    tfParams.earQ = 8;
    tfParams.agc = 1;
    tfParams.differ = 1;
    tfParams.tau = 3;
    tfParams.step = 0.5;
    
    preprocStimParams.tfParams = tfParams;
    preprocStimParams.outputDir = preprocDir;
    lpat = sprintf('earQ_%d.agc_%d.differ_%d', tfParams.earQ, tfParams.agc, tfParams.differ);
    preprocStimParams.outputPattern = ['stim.preproc-' preprocType '-' lpat '.%d.mat'];
    
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
    stimInfo.zscored = zscore;
    respInfo.zscored = zscore;
    if zscore
        fprintf('z-scoring stimulus and response...\n');
        wholeStim = norm_std_mean(wholeStim);   
        respInfo.mean = mean(wholeResponse);
        respInfo.std = std(wholeResponse);
        wholeResponse = (wholeResponse - respInfo.mean) / respInfo.std;
    end
    %}
    
    
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
    