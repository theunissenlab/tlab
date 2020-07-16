function srdata = get_sr_data(dataRootDir, stimsDir, cellName)

    cellDir = fullfile(dataRootDir, cellName);

    %% get datafiles in cell directory
    datasets = find_datasets(cellDir, stimsDir, 'conspecific|zfsongs');    
    stimFiles = datasets{1}.srPairs.stimFiles;
    respFiles = datasets{1}.srPairs.respFiles;
    dataDir = datasets{1}.dirname;
    
    %% preprocess the sound to produce spectrograms
    preprocType = 'stft';
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
    preprocStimParams.outputPattern = ['stim.preproc-' preprocType '.%d.mat'];
    
    [wholeStim, groupIndex, stimInfo, preprocStimParams] = preprocSound(stimFiles, preprocStimParams);

    %% z-score the stimulus
    zscore = 1;
    if zscore
        fprintf('z-scoring stimulus...\n');
        wholeStim = norm_std_mean(wholeStim);    
    end
    
    
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

    
    %% split response
    preprocRespParams.split = 1;
    [wholeSplitResponse, respSplitInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);
    
    srdata.stim = wholeStim;
    srdata.groupIndex = groupIndex;
    srdata.stimInfo = stimInfo;
    srdata.stimParams = preprocStimParams;
    srdata.resp = wholeResponse;
    srdata.respInfo = respInfo;
    srdata.respParams = preprocRespParams;
    srdata.splitResp = wholeSplitResponse;
    

    