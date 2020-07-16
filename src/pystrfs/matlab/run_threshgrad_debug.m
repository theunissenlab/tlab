function run_threshgrad_debug(unitFile, preprocFile, stimClass, trainingGroups, validationGroups, strfLength, transforms, threshold, outputFile)

    if exist(outputFile, 'file')
        fprintf('Output file already exists, returning...\n');
        return;
    end

    stimRespData = get_stimresp_data(unitFile, preprocFile, stimClass);
    
    %% transform stimulus
    wholeStim = stimRespData.wholeStim;    
    for k = 1:length(transforms)       
        fprintf('Performing %s transform...\n', transforms{k});
        tfunc = sprintf('transform_%s(wholeStim)', transforms{k});
        wholeStim = eval(tfunc);        
    end    
    stimRespData.wholeStim = wholeStim;
       
    %% Initialize strflab global variables with our stim and responses
    global globDat;
    strfData(stimRespData.wholeStim, stimRespData.wholeResp, stimRespData.groupIndex);

    
    %% Initialize a linear model
    strfDelays = 0:(strfLength-1);
    modelParams = linInit(stimRespData.numChannels, strfDelays, 'linear');
    
    
    %% Initialize Threshold Gradient Descent options
    optOptions = trnThreshGradDescLS();        
    optOptions.threshold = threshold;    
    optOptions.earlyStop = 0;
    optOptions.display = 1;
    optOptions.debug = 1;
    
    optOptions.maxIter = 1000;
    optOptions.maxStep = 1;
    optOptions.maxLSIter = 50;
    optOptions.stepSize = 1e-4;
    optOptions.stepConvWindow = 8;
    optOptions.stepConv = 1e-10;
    
    %% compute STRF across entire dataset without early stopping
    trainingIndex = findIdx(trainingGroups, stimRespData.groupIndex);
    validationIndex = findIdx(validationGroups, stimRespData.groupIndex);    
    
    [modelParamsTrained, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);
    
    save(outputFile, 'unitFile', 'preprocFile', 'stimClass', 'trainingGroups', 'validationGroups', 'strfLength', 'transforms', ...
                     'threshold', 'stimRespData', 'strfDelays', 'modelParams', 'optOptions', 'trainingIndex', 'validationIndex', ...
                     'modelParamsTrained');