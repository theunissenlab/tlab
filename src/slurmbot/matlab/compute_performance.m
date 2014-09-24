function perfData = compute_performance(preprocData, trainingIndex, validationIndex, modelParamsTrained)

    global globDat;

    stimInfo = preprocData.stimInfo;
    respInfo = preprocData.respInfo;
    wholeSplitResponse = preprocData.splitResp;
    avgNumTrials = mean(preprocData.respInfo.numTrials); %taking mean isn't necessary here
    infoFreqCutoff = 90; %Hz
    infoWindowSize = 0.250; %250ms
    
    
    %% un-zscore real response
    wholeResp = globDat.resp;
    if respInfo.zscored
        wholeResp = wholeResp*respInfo.std + respInfo.mean;
    end
    
    %% Get concatenated responses of training and early stopping responses 
    respTrainReal = wholeResp(trainingIndex);
    respTrainRealHalf1 = wholeSplitResponse(1, trainingIndex);
    respTrainRealHalf2 = wholeSplitResponse(2, trainingIndex);
    
    respValidReal = wholeResp(validationIndex);
    respValidRealHalf1 = wholeSplitResponse(1, validationIndex);
    respValidRealHalf2 = wholeSplitResponse(2, validationIndex);
    
    %% Compute model responses, coherences, and un-zscore them
    if ~isempty(trainingIndex)
        [modelParamsTrained, respTrain] = strfFwd(modelParamsTrained, trainingIndex);
        respTrain(isnan(respTrain)) = 0;
        if respInfo.zscored
            respTrain = respTrain*respInfo.std + respInfo.mean;
        end
        respTrain(respTrain < 0) = 0;
        respTrain = (respTrain / max(respTrain))*max(respTrainReal);
    
        [cBoundTrain, cTrain] = compute_coherence_full(respTrain, respTrainReal, respTrainRealHalf1, respTrainRealHalf2, stimInfo.sampleRate, avgNumTrials, infoFreqCutoff, infoWindowSize);
        trainingPerfRatio = cTrain.info / cBoundTrain.info;    
    else
        respTrain = -1;
        cBoundTrain = struct;
        cTrain = struct;
        trainingPerfRatio = -1;
    end
    
    if ~isempty(validationIndex)    
        [modelParamsTrained, respValid] = strfFwd(modelParamsTrained, validationIndex);    
        respValid(isnan(respValid)) = 0;
        if respInfo.zscored
            respValid = respValid*respInfo.std + respInfo.mean;
        end
        respValid(respValid < 0) = 0;
        respValid = (respValid / max(respValid))*max(respValidReal);
    
        [cBoundValid, cValid] = compute_coherence_full(respValid, respValidReal, respValidRealHalf1, respValidRealHalf2, stimInfo.sampleRate, avgNumTrials, infoFreqCutoff, infoWindowSize);
        validPerfRatio = cValid.info / cBoundValid.info;
    else       
        respValid = -1;
        cBoundValid = struct;
        cValid = struct;
        validPerfRatio = -1;
    end
    
    perfData = struct;
    perfData.wholeResp = wholeResp;
    perfData.train.resp = respTrain;
    perfData.train.cBound = cBoundTrain;
    perfData.train.c = cTrain;
    perfData.train.perf = trainingPerfRatio;
    
    perfData.valid.resp = respValid;
    perfData.valid.cBound = cBoundValid;
    perfData.valid.c = cValid;
    perfData.valid.perf = validPerfRatio;
    