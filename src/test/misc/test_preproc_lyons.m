function test_preproc_lyons()

    %cellName = 'pupu2122_2_B';
    %cellName = 'yy1617_4_A';
    %cellName = 'gg0116_1_A';
    %cellName = 'obla1305_7';
    %cellName = 'yy1617_5_B';
    cellName = 'oo0909_2_A';

    dataRootDir = '/auto/fdata/mschachter/data';
    stimsDir = fullfile(dataRootDir, 'all_stims');
    cellDir = fullfile(dataRootDir, cellName);
    
    trainingGroups = 1:18;
    validationGroups = 19:20;
    
    strfLength = 75;
    nlType = 'linear';

    %% preproc stim and responses
    preprocData = preproc_stim_response_lyons(cellDir, stimsDir, 0);    
    groupIndex = preprocData.groupIndex;
    wholeStim = preprocData.stim;
    wholeResponse = preprocData.resp;
    dataDir = preprocData.dataDir;
    stimInfo = preprocData.stimInfo;
    
    displayStimResp = 0;
    if displayStimResp       
        display_stim_resp(preprocData);
    end
    
    
    %% get training and early stopping dataset indicies
    trainingIndex = findIdx(trainingGroups, groupIndex);
    validationIndex = findIdx(validationGroups, groupIndex);

    %% Initialize strflab global variables with our stim and responses
    global globDat
    strfData(wholeStim, wholeResponse, groupIndex);


    %% Initialize a linear model
    strfDelays = 0:(strfLength-1);
    modelParams = linInit(stimInfo.numStimFeatures, strfDelays, nlType);
    modelParams.b1 = mean(wholeResponse);
    
    
    %% Initialize Direct Fit options
    optOptions = trnDirectFit();
    optOptions.tolerances = [0.500 0.100 0.050 0.010 0.005 1e-03 5e-04 1e-04 5e-05];
    optOptions.sparsenesses = [0 1 2 4 6 8 10 12];
    optOptions.infoFreqCutoff = 90;
    optOptions.infoWindowSize = 0.250;
       
    %% run direct fit
    [modelParamsTrained, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);
    
    %% reset the response data to it's non-inverted form
    strfData(wholeStim, wholeResponse, groupIndex);
    
    perfData = compute_performance(preprocData, trainingIndex, validationIndex, modelParamsTrained);
                     
    %% display stuff
    display = 1;
    if display
        display_perfdata(preprocData, perfData, modelParamsTrained, trainingIndex, validationIndex);
    end
    %}