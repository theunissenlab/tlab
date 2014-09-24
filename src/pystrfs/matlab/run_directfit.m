function run_directfit(unitFile, preprocFile, stimClass, trainingGroups, validationGroups, strfLength, transforms, outputFile)

    stimRespData = get_stimresp_data(unitFile, preprocFile, stimClass);
    
    %% transform stimulus
    wholeStim = stimRespData.wholeStim;
    
    for k = 1:length(transforms)       
        fprintf('Performing %s transform...\n', transforms{k});
        tfunc = sprintf('transform_%s(wholeStim)', transforms{k});
        wholeStim = eval(tfunc);        
    end
    
    stimRespData.wholeStim = wholeStim;
    
    
    %% get training and early stopping dataset indicies
    trainingIndex = findIdx(trainingGroups, stimRespData.groupIndex);
    validationIndex = findIdx(validationGroups, stimRespData.groupIndex);
    
    
    %% Initialize strflab global variables with our stim and responses
    global globDat;
    strfData(stimRespData.wholeStim, stimRespData.wholeResp, stimRespData.groupIndex);

    
    %% Initialize a linear model
    strfDelays = 0:(strfLength-1);
    modelParams = linInit(stimRespData.numChannels, strfDelays, 'linear');
    
    
    %% Initialize Direct Fit options
    optOptions = trnDirectFit();
    optOptions.tolerances = [0.100 0.050 0.010 0.005 1e-03 5e-04 1e-04];
    optOptions.sparsenesses = [0 1 2 3 4 5 6 7];
    optOptions.infoFreqCutoff = 90;
    optOptions.infoWindowSize = 0.250;
    optOptions.biasForHighTol = 0;
    
    
    %% run direct fit
    [modelParamsTrained, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);
    
    
    %% produce responses
    [modelResp, rawModelResp] = compute_response(stimRespData, modelParamsTrained);
    
    if exist(outputFile, 'file')
        delete(outputFile);
    end
    
    trStrs = '';
    for k = 1:length(transforms)
        if k > 1
            ts = ',';
        else
            ts = '';
        end
        trStrs = [trStrs ts transforms{k}];
    end
    
    if isempty(trStrs)
        trStrs = ' ';
    end
    
    h5 = h5utils();
    fid = h5.create(outputFile);    
    
    h5.set_attr(fid, '/data', 'sample_rate', stimRespData.sampleRate);
    h5.set_attr(fid, '/data', 'unit_file', unitFile);
    h5.set_attr(fid, '/data', 'preproc_file', preprocFile);
    h5.set_attr(fid, '/data', 'stim_class', stimClass);    
    
    h5.set_ds(fid, '/data', 'training_groups', trainingGroups);
    h5.set_ds(fid, '/data', 'validation_groups', validationGroups);
    h5.set_ds(fid, '/data', 'group_index', stimRespData.groupIndex);
    h5.set_ds(fid, '/data', 'training_index', trainingIndex);
    h5.set_ds(fid, '/data', 'validation_index', validationIndex);
    
    h5.set_attr(fid, '/opt', 'method', 'direct_fit');
    h5.set_ds(fid, '/opt', 'tolerances', optOptions.tolerances);
    h5.set_ds(fid, '/opt', 'sparsenesses', optOptions.sparsenesses);
    h5.set_ds(fid, '/opt', 'info_freq_cutoff', optOptions.infoFreqCutoff);
    h5.set_ds(fid, '/opt', 'info_window_size', optOptions.infoWindowSize);
        
    h5.set_attr(fid, '/model', 'type', 'linear');
    h5.set_attr(fid, '/model', 'strf_length', strfLength);
    h5.set_attr(fid, '/model', 'transforms', trStrs);
    h5.set_ds(fid, '/model', 'weights', modelParamsTrained.w1);
    h5.set_ds(fid, '/model', 'bias', modelParamsTrained.b1);
    h5.set_ds(fid, '/model', 'response', modelResp);
    h5.set_ds(fid, '/model', 'response_raw', rawModelResp);
    
    h5.close(fid);
    