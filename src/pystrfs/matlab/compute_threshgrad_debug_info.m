function compute_threshgrad_debug_info(debugFile)

    vars = load(debugFile);

    stimRespData = vars.stimRespData;
    wts = vars.optOptions.diagnostics.params;
    grads = vars.optOptions.diagnostics.grads;
    modelParamsTrained = vars.modelParamsTrained;        
    validationIndex = vars.validationIndex;    
    thresh = vars.threshold;
    clear vars;
    
    wholeResp = stimRespData.wholeResp;
    wholeSplitResp = stimRespData.wholeSplitResp;
    numTrials = stimRespData.numTrials;
    sampleRate = stimRespData.sampleRate;
    
    infos = zeros(length(wts), 1);
    errs = zeros(length(wts), 1);
    
    gradMax = zeros(length(wts), 1);
    gradMean = zeros(length(wts), 1);
    gradStd = zeros(length(wts), 1);    
    gradMaxThresh = zeros(length(wts), 1);
    gradMeanThresh = zeros(length(wts), 1);
    gradStdThresh = zeros(length(wts), 1);    
    
    global globDat;
    strfData(stimRespData.wholeStim, stimRespData.wholeResp, stimRespData.groupIndex);
    
    %compute validation performance for each STRF
    for k = 1:length(wts)    
        
        %% compute grad stats
        g = grads{k};
        ag = abs(g);
        [mg, ig] = max(ag);
                
        tg = g(ag >= thresh*mg);
        atg = abs(tg);
        
        gradMax(k) = g(ig);
        gradMean(k) = mean(ag);
        gradStd(k) = std(ag);
        
        [amg, aig] = max(atg);
        gradMaxThresh(k) = tg(aig);
        gradMeanThresh(k) = mean(atg);
        gradStdThresh(k) = std(atg);
        
        %% compute response
        mp = strfUnpak(modelParamsTrained, wts{k});
        [modelResp, ~] = compute_response(stimRespData, mp);
        
        %% Get concatenated responses of validation responses        
        respValid = modelResp(validationIndex);
        respValidReal = wholeResp(validationIndex);
        respValidRealHalf1 = wholeSplitResp(1, validationIndex);
        respValidRealHalf2 = wholeSplitResp(2, validationIndex);

        avgNumTrials = mean(numTrials);

        infoFreqCutoff = 90;
        infoWindowSize = 0.200;

        %% Compute coherences and msq error
        [cBoundValid, cValid] = compute_coherence_full(respValid, respValidReal, respValidRealHalf1, respValidRealHalf2, sampleRate, avgNumTrials, infoFreqCutoff, infoWindowSize);
        infos(k) = cValid.info;
        errs(k) = sum((respValid' - wholeResp(validationIndex)) .^ 2);
        
        fprintf('Iteration %d: err=%0.3f, info=%0.2f\n', k, errs(k), infos(k));
    end
    
    infoDebugFile = sprintf('%s.info.mat', debugFile);
    fprintf('Saving info file at %s\n', infoDebugFile);
    save(infoDebugFile, 'infos', 'errs', 'gradMax', 'gradMean', 'gradStd', 'gradMaxThresh', 'gradMeanThresh', 'gradStdThresh');
    