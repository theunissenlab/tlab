function compute_response_info_coherence(responseFile)

    %% get parameters from optimization output file
    h5 = h5utils();
    fid = h5.open(responseFile, 'a');
    
    unitFileName = h5.get_attr(fid, '/data', 'unit_file');
    unit = read_unit_file(unitFileName);
    
    stimClass = h5.get_attr(fid, '/data', 'stim_class');
    sampleRate = h5.get_attr(fid, '/data', 'sample_rate');
    
    groupIndex = h5.get_ds(fid, '/data/group_index');
    trainingIndex = h5.get_ds(fid, '/data/training_index');
    validationIndex = h5.get_ds(fid, '/data/validation_index');
        
    modelResp = h5.get_ds(fid, '/model/response');
    
    %% get spike data from unit
    allResps = unit.class_responses.(stimClass);    
    try
        stimIndices = h5.get_ds(fid, '/data/stim_indices_used');
    catch
        stimIndices = 1:len(allResps);
    end
    
    allSpikeTrials = cell(length(stimIndices), 1);
    numTrials = zeros(length(stimIndices), 1);
    respLengths = zeros(length(stimIndices), 1);
    rindex = 1;
    for si = 1:length(stimIndices)
        k = stimIndices(si);        
        allSpikeTrials{rindex} = allResps{k}.spikeTrials;        
        numTrials(rindex) = allResps{k}.numTrials;
        respLengths(rindex) = sum(groupIndex == rindex) / sampleRate;
        rindex = rindex + 1;
    end
    
    
    %% compute PSTH and split PSTHs
    preprocRespParams = struct;         %create preprocessing struct
    preprocRespParams.units = 's';      %files have units of s
    preprocRespParams.binSize = 0.001;  %1ms bin size, same sample rate as spectrograms
    preprocRespParams.stimLengths = respLengths; %needed to compute correct PSTH lengths

    [wholeResp, groupIndex, respInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);

    preprocRespParams.split = 1;
    [wholeSplitResp, respSplitInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);

    
    %% get model responses    
    respTrain = modelResp(trainingIndex);
    respValid = modelResp(validationIndex);
    
    
    %% Get concatenated responses of training and early stopping responses 
    respTrainReal = wholeResp(trainingIndex);
    respTrainRealHalf1 = wholeSplitResp(1, trainingIndex);
    respTrainRealHalf2 = wholeSplitResp(2, trainingIndex);
    
    respValidReal = wholeResp(validationIndex);
    respValidRealHalf1 = wholeSplitResp(1, validationIndex);
    respValidRealHalf2 = wholeSplitResp(2, validationIndex);
    
    avgNumTrials = mean(numTrials);
    
    infoFreqCutoff = 90;
    infoWindowSize = 0.200;
    
    %% Compute model responses, coherences, and un-zscore them
    if ~isempty(trainingIndex)        
        [cBoundTrain, cTrain] = compute_coherence_full(respTrain, respTrainReal, respTrainRealHalf1, respTrainRealHalf2, sampleRate, avgNumTrials, infoFreqCutoff, infoWindowSize);
        trainingPerfRatio = cTrain.info / cBoundTrain.info;    
    else        
        cBoundTrain = struct;
        cTrain = struct;
        trainingPerfRatio = -1;
    end
    
    if ~isempty(validationIndex)    
        [cBoundValid, cValid] = compute_coherence_full(respValid, respValidReal, respValidRealHalf1, respValidRealHalf2, sampleRate, avgNumTrials, infoFreqCutoff, infoWindowSize);
        validPerfRatio = cValid.info / cBoundValid.info;
    else       
        cBoundValid = struct;
        cValid = struct;
        validPerfRatio = -1;
    end
    
    
    %% print summary data
    fprintf('Coherence performance info:\n');
    fprintf('\tTraining Upper Bound: %0.1f; CI=(%0.2f, %0.2f) bits\n', cBoundTrain.info, cBoundTrain.infoLower, cBoundTrain.infoUpper);
    fprintf('\tTraining: %0.f bits, ratio=%0.3f\n', cTrain.info, trainingPerfRatio);
    fprintf('\tValidation Upper Bound: %0.1f; CI=(%0.2f, %0.2f) bits\n', cBoundValid.info, cBoundValid.infoLower, cBoundValid.infoUpper);
    fprintf('\tValidation: %0.f bits, ratio=%0.3f\n', cValid.info, validPerfRatio);
    
    
    %% write to file
    h5.set_ds(fid, '/model/performance/coherence/bound/training', 'frequency', cBoundTrain.f);
    h5.set_ds(fid, '/model/performance/coherence/bound/training', 'coherence_mean', cBoundTrain.c);
    h5.set_ds(fid, '/model/performance/coherence/bound/training', 'coherence_upper', cBoundTrain.cUpper);
    h5.set_ds(fid, '/model/performance/coherence/bound/training', 'coherence_lower', cBoundTrain.cLower);
    h5.set_ds(fid, '/model/performance/coherence/bound/training', 'info_mean', cBoundTrain.info);
    h5.set_ds(fid, '/model/performance/coherence/bound/training', 'info_upper', cBoundTrain.infoUpper);
    h5.set_ds(fid, '/model/performance/coherence/bound/training', 'info_lower', cBoundTrain.infoLower);
    
    h5.set_ds(fid, '/model/performance/coherence/bound/validation', 'frequency', cBoundValid.f);
    h5.set_ds(fid, '/model/performance/coherence/bound/validation', 'coherence_mean', cBoundValid.c);
    h5.set_ds(fid, '/model/performance/coherence/bound/validation', 'coherence_upper', cBoundValid.cUpper);
    h5.set_ds(fid, '/model/performance/coherence/bound/validation', 'coherence_lower', cBoundValid.cLower);
    h5.set_ds(fid, '/model/performance/coherence/bound/validation', 'info_mean', cBoundValid.info);
    h5.set_ds(fid, '/model/performance/coherence/bound/validation', 'info_upper', cBoundValid.infoUpper);
    h5.set_ds(fid, '/model/performance/coherence/bound/validation', 'info_lower', cBoundValid.infoLower);
    
    h5.set_ds(fid, '/model/performance/coherence/training', 'frequency', cTrain.f);
    h5.set_ds(fid, '/model/performance/coherence/training', 'coherence_mean', cTrain.c);
    h5.set_ds(fid, '/model/performance/coherence/training', 'coherence_upper', cTrain.cUpper);
    h5.set_ds(fid, '/model/performance/coherence/training', 'coherence_lower', cTrain.cLower);
    h5.set_ds(fid, '/model/performance/coherence/training', 'info_mean', cTrain.info);
    h5.set_ds(fid, '/model/performance/coherence/training', 'info_upper', cTrain.infoUpper);
    h5.set_ds(fid, '/model/performance/coherence/training', 'info_lower', cTrain.infoLower);
    
    h5.set_ds(fid, '/model/performance/coherence/validation', 'frequency', cValid.f);
    h5.set_ds(fid, '/model/performance/coherence/validation', 'coherence_mean', cValid.c);
    h5.set_ds(fid, '/model/performance/coherence/validation', 'coherence_upper', cValid.cUpper);
    h5.set_ds(fid, '/model/performance/coherence/validation', 'coherence_lower', cValid.cLower);
    h5.set_ds(fid, '/model/performance/coherence/validation', 'info_mean', cValid.info);
    h5.set_ds(fid, '/model/performance/coherence/validation', 'info_upper', cValid.infoUpper);
    h5.set_ds(fid, '/model/performance/coherence/validation', 'info_lower', cValid.infoLower);
    
    
    h5.close(fid);
    