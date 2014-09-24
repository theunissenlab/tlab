function test_preproc_whiten()

    cellName = 'pupu2122_2_B';
    %cellName = 'yy1617_4_A';

    dataRootDir = '/auto/fdata/mschachter/data';
    stimsDir = fullfile(dataRootDir, 'all_stims');
    cellDir = fullfile(dataRootDir, cellName);
    
    trainingGroups = 1:18;
    validationGroups = 19:20;
    
    strfLength = 75;
    nlType = 'linear';

    %% preproc stim and responses
    preprocData = preproc_stim_response(cellDir, stimsDir, 0);    
    groupIndex = preprocData.groupIndex;
    wholeStim = preprocData.stim;
    wholeResponse = preprocData.resp;
    dataDir = preprocData.dataDir;
    stimInfo = preprocData.stimInfo;
    
    %% whiten stims
    wparams = struct;
    wparams.method = 'meansub';
    [whitenedStim, wparams, winfo] = preprocWhiten(wholeStim, wparams);    
    
    %{
    uindx = unique(groupIndex);
    for k = uindx
       
        g = uindx(k);
        gindx = find(groupIndex == g);
        stim = wholeStim(gindx, :);
        wstim = whitenedStim(gindx, :);
        
        figure; hold on;
        subplot(2, 1, 1); hold on;
        imagesc(stim'); axis tight;
        colorbar;
        subplot(2, 1, 2); hold on;
        imagesc(wstim'); axis tight;
        colorbar;
        
    end
    %}
    
    wholeStim = whitenedStim;
    
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
    