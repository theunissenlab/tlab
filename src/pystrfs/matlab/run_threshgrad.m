function run_threshgrad(unitFile, preprocFile, stimClass, trainingGroups, validationGroups, strfLength, transforms, threshold, outputFile, outputNL)

    if nargin < 10
        outputNL = 'linear';
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
    
    %% transform response
    if strcmp(outputNL, 'exponential')       
        wholeResp = stimRespData.wholeResp;
        gindx = stimRespData.groupIndex;
        grps = unique(gindx);
        for k = 1:length(grps)           
            gi = gindx == grps(k);
            wholeResp(gi) = wholeResp(gi)*stimRespData.numTrials(k);            
        end
        stimRespData.wholeResp = wholeResp;
    end
    
    
    %% get training and early stopping dataset indicies
    trainingIndex = findIdx(trainingGroups, stimRespData.groupIndex);
    validationIndex = findIdx(validationGroups, stimRespData.groupIndex);
        
    
    %% Initialize strflab global variables with our stim and responses
    global globDat;
    strfData(stimRespData.wholeStim, stimRespData.wholeResp, stimRespData.groupIndex);

    
    %% Initialize a linear model
    strfDelays = 0:(strfLength-1);
    modelParams = linInit(stimRespData.numChannels, strfDelays, outputNL);
    
    
    %% Initialize Direct Fit options
    %optOptions = trnThreshGradDescLS();    
    optOptions = trnThreshGradDesc();    
    optOptions.lineSearch = '';
    optOptions.stepSize = 1e-4;
    optOptions.threshold = threshold;
    optOptions.earlyStop = 0;
    optOptions.display = 1;
    
    %% run direct fit
    [modelParamsTrained, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);
    
    
    %% produce responses
    [modelResp, rawModelResp] = compute_response(stimRespData, modelParamsTrained);
    if strcmp(outputNL, 'exponential')       
        modelResp = rawModelResp;
    end
    
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
    
    h5.set_ds(fid, '/data', 'group_index', stimRespData.groupIndex);
    h5.set_ds(fid, '/data', 'training_index', trainingIndex);
    h5.set_ds(fid, '/data', 'validation_index', validationIndex);
    
    h5.set_attr(fid, '/opt', 'method', 'threshgrad');
    h5.set_ds(fid, '/opt', 'tolerances', optOptions.tolerances);
    h5.set_ds(fid, '/opt', 'sparsenesses', optOptions.sparsenesses);
    h5.set_ds(fid, '/opt', 'info_freq_cutoff', optOptions.infoFreqCutoff);
    h5.set_ds(fid, '/opt', 'info_window_size', optOptions.infoWindowSize);
        
    h5.set_attr(fid, '/model', 'strf_length', strfLength);
    h5.set_attr(fid, '/model', 'transforms', trStrs);
    h5.set_attr(fid, '/model', 'output_nl', outputNL);
    h5.set_ds(fid, '/model', 'weights', modelParamsTrained.w1);
    h5.set_ds(fid, '/model', 'bias', modelParamsTrained.b1);
    h5.set_ds(fid, '/model', 'response', modelResp);
    h5.set_ds(fid, '/model', 'response_raw', rawModelResp);
    
    h5.close(fid);
    