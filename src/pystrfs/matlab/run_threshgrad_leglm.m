function run_threshgrad_leglm(unitFile, preprocFile, stimClass, trainingGroups, validationGroups, strfLength, transforms, threshold, outputFile, familyName, outputNL)

    stimRespData = get_stimresp_data(unitFile, preprocFile, stimClass);
    
    numTrials = round(mean(stimRespData.numTrials));
    
    %% set up GLM family
    switch familyName
        case 'gaussian'
            family = init_glm_family_gaussian();
        case 'binomial'
            family = init_glm_family_binomial(numTrials);
        case 'poisson'
            family = init_glm_family_poisson();
    end    
    
    %% transform stimulus
    wholeStim = stimRespData.wholeStim;    
    for k = 1:length(transforms)       
        fprintf('Performing %s transform...\n', transforms{k});
        tfunc = sprintf('transform_%s(wholeStim)', transforms{k});
        wholeStim = eval(tfunc);        
    end    
    stimRespData.wholeStim = wholeStim;
        
    %% transform response
    if strcmp(family.type, 'poisson')
        wholeResp = stimRespData.wholeResp;
        gindx = stimRespData.groupIndex;
        grps = unique(gindx);
        for k = 1:length(grps)           
            gi = gindx == grps(k);
            wholeResp(gi) = wholeResp(gi)*stimRespData.numTrials(k);            
        end
        stimRespData.wholeResp = wholeResp;
    end
    
    %% Initialize strflab global variables with our stim and responses
    global globDat;
    strfData(stimRespData.wholeStim, stimRespData.wholeResp, stimRespData.groupIndex);

    
    %% Initialize a linear model
    strfDelays = 0:(strfLength-1);
    dispersion = 1;
    modelParams = leglmInit(stimRespData.numChannels, strfDelays, family, dispersion);
    
    
    %% Initialize Threshold Gradient Descent options
    optOptions = trnThreshGradDescLS();
    optOptions.threshold = threshold;
    optOptions.earlyStop = 1;
    optOptions.display = 1;   
    
    optOptions.maxLSIter = 100;
    optOptions.maxIter = 500;
    optOptions.errLastN = 20;
    optOptions.errStartN = 2;
    optOptions.stepConvWindow = 8;
    optOptions.maxStep = 5;
    optOptions.minStep = 1e-6;
    optOptions.stepSize = 1e-2;
    optOptions.stepConv = 1e-5;
            
    if length(trainingGroups) == 18
        nPrts = 3;
    elseif length(trainingGroups) == 9
        nPrts = 3;
    end
    
    prts = cv_partition(nPrts, trainingGroups);
    partitions = {};
    for k = 1:length(prts)
        if ~isempty(prts{k}.validation)
            partitions{end+1} = prts{k};
        end
    end
    
    strfs = zeros(length(partitions), stimRespData.numChannels, strfLength);
    biases = zeros(length(partitions), 1);
    numIters = zeros(length(partitions), 1);
    runTimes = zeros(length(partitions), 1);
    holdOuts = zeros(length(partitions), round(length(trainingGroups) / nPrts));
    
    for k = 1:length(partitions)
       
        tic;
        tg = partitions{k}.training;
        eg = partitions{k}.validation;
        holdOuts(k, :) = eg;
        
        %% get training and early stopping dataset indicies
        trIndex = findIdx(tg, stimRespData.groupIndex);
        esIndex = findIdx(eg, stimRespData.groupIndex);        
        
        %% run direct fit
        [modelParamsTrained, optOptions] = strfOpt(modelParams, trIndex, optOptions, esIndex);
        etime = toc;
        fprintf('\nCV Run #%d took %f\n', k, etime);
        
        strfs(k, :, :) = squeeze(modelParamsTrained.w1);
        biases(k) = modelParamsTrained.b1;
        numIters(k) = optOptions.diagnostics.bestiter;
        runTimes(k) = etime;
    end
    
    numIters
    iterMean = mean(numIters);
    iterStd = std(numIters);
    fprintf('# of iterations: %0.0f +/- %0.0f\n', iterMean, iterStd);
    
    %% compute STRF across entire dataset without early stopping
    trainingIndex = findIdx(trainingGroups, stimRespData.groupIndex);
    validationIndex = findIdx(validationGroups, stimRespData.groupIndex);
    optOptions.earlyStop = 0;
    optOptions.maxIter = iterMean;
    [modelParamsTrained, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);
    
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
    
    %% produce responses
    [modelResp, rawModelResp] = compute_response(stimRespData, modelParamsTrained, 0);
    
    h5 = h5utils();
    fid = h5.create(outputFile);    
    
    h5.set_attr(fid, '/data', 'sample_rate', stimRespData.sampleRate);
    h5.set_attr(fid, '/data', 'unit_file', unitFile);
    h5.set_attr(fid, '/data', 'preproc_file', preprocFile);
    h5.set_attr(fid, '/data', 'stim_class', stimClass);    
    
    h5.set_ds(fid, '/data', 'group_index', stimRespData.groupIndex);
    h5.set_ds(fid, '/data', 'training_index', trainingIndex);
    h5.set_ds(fid, '/data', 'validation_index', validationIndex);
    h5.set_ds(fid, '/data', 'training_groups', trainingGroups);
    h5.set_ds(fid, '/data', 'validation_groups', validationGroups);
    
    h5.set_ds(fid, '/data', 'holdout_subgroups', holdOuts);
    
    h5.set_attr(fid, '/opt', 'method', 'threshgrad');
    h5.set_ds(fid, '/opt', 'threshold', optOptions.threshold);
    h5.set_ds(fid, '/opt', 'max_iter', optOptions.maxIter);
    h5.set_ds(fid, '/opt', 'num_iters', numIters);
    h5.set_ds(fid, '/opt', 'run_times', runTimes);
            
    h5.set_attr(fid, '/model', 'family_name', familyName);
    h5.set_attr(fid, '/model', 'strf_length', strfLength);
    h5.set_attr(fid, '/model', 'transforms', trStrs);
    h5.set_ds(fid, '/model', 'cv_weights', strfs);
    h5.set_ds(fid, '/model', 'cv_bias', biases);
    h5.set_ds(fid, '/model', 'weights', squeeze(modelParamsTrained.w1));
    h5.set_ds(fid, '/model', 'bias', modelParamsTrained.b1);    
    h5.set_ds(fid, '/model', 'm', modelParamsTrained.m);    
    h5.set_ds(fid, '/model', 'response', modelResp);
    h5.set_ds(fid, '/model', 'response_raw', rawModelResp);
        
    h5.close(fid);
end

function [y, dy] = exp_outputnl(x)
    y = exp(x);
    if nargout > 1
        dy = exp(x);
    end
end

function [y, dy] = sigmoid_outputnl(x)
    y = 1 ./ (1 + exp(-x));
    if nargout > 1
        dy = y .* (1 - y);
    end
end

function [y, dy] = identity_outputnl(x)
    y = x;
    if nargout > 1
        dy = 1;
    end
end