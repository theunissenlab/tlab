function run_threshgrad_cv(unitFile, preprocFile, stimClass, trainingGroups, validationGroups, earlyStoppingGroups, strfLength, transforms, threshold, outputFile, modelType)

    debug = 0;

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
    if ismember(modelType, {'poisson', 'leglm'})
        fprintf('Transforming response to spike counts...\n');
        wholeResp = stimRespData.wholeResp;
        gindx = stimRespData.groupIndex;
        grps = unique(gindx);
        for k = 1:length(grps)           
            gi = gindx == grps(k);
            wholeResp(gi) = wholeResp(gi)*stimRespData.numTrials(k);            
        end
        stimRespData.wholeResp = wholeResp;
    end
    
    %display_stim_resp_data(stimRespData);
    
    %% Initialize strflab global variables with our stim and responses
    global globDat;
    strfData(stimRespData.wholeStim, stimRespData.wholeResp, stimRespData.groupIndex);

    
    %% Initialize a linear model
    strfDelays = 0:(strfLength-1);
    if strcmp('linear', modelType)
        modelParams = linInit(stimRespData.numChannels, strfDelays, modelType);
    elseif ismember(modelType, {'poisson', 'binomial', 'leglm'})
        modelParams = init_glm_model(stimRespData, strfDelays, modelType);
    else
        error('Unknown model type %s', modelType);
    end
    
    
    %% Initialize Threshold Gradient Descent options
    optOptions = trnThreshGradDescLS();        
    optOptions.earlyStop = 1;
    optOptions.lineSearch = 1;
    optOptions.display = 1;
    optOptions.debug = debug;
    
    optOptions.maxIter = 2500;    
    
    optOptions.stepSize = 1e-4;    
    optOptions.stepConv = 0;
    if threshold < 1
        optOptions.maxStep = 1;
        optOptions.maxLSIter = 50;
    else
        optOptions.maxStep = 10;
        optOptions.maxLSIter = 75;
    end

    switch threshold
        case 0.0
            optOptions.stepSize = 0.0006;
        case 0.25
            optOptions.stepSize = 0.0006;
        case 0.50
            optOptions.stepSize = 0.00065;
        case 0.75
            optOptions.stepSize = 0.00075;
        case 1.00
            optOptions.stepSize = 0.15;
    end
    %fprintf('Using step size %f\n', optOptions.stepSize);
    
    %% set up param groups, ensure that bias is not thresholded out
    if isempty(strfind(preprocFile, 'surprise'))               
        optOptions.paramGroups = ones(1, stimRespData.numChannels*length(strfDelays) + 1);
        optOptions.paramGroups(end) = 2;        
        optOptions.threshold = [threshold 0];    
    else
        %% split the threshold into three parts if we're doing surprise
        nhalf = round(stimRespData.numChannels / 2);
        nhalfParams = nhalf * length(strfDelays);
        paramGroups = ones(1, 2*nhalfParams + 1);
        paramGroups(nhalfParams+1:(end-1)) = 2;
        paramGroups(end) = 3;
        
        optOptions.paramGroups = paramGroups;
        optOptions.threshold = [threshold threshold 0];        
    end
    
    if strcmp(modelType, 'leglm')       
        %% add extra parameter for m
        optOptions.paramGroups(end+1) = max(optOptions.paramGroups) + 1;        
        optOptions.threshold(end+1) = 0;
    end
    
    %% split the early stopping groups
    numEsGroups = 3;
    esLen = length(earlyStoppingGroups);
    esPerTrain = round(esLen / numEsGroups);
    esGroups = earlyStoppingGroups(reshape(1:esLen, esPerTrain, numEsGroups)');
    
    strfs = zeros(numEsGroups, stimRespData.numChannels, strfLength);
    biases = zeros(numEsGroups, 1);
    numIters = zeros(numEsGroups, 1);
    runTimes = zeros(numEsGroups, 1);
    if strcmp(modelType, 'leglm')
        ms = zeros(numEsGroups, 1);
    end
    
    totalTime = 0;
    
    %% train for each early stopping group and produce a STRF
    for k = 1:numEsGroups 
        eg = esGroups(k, :);
        tg = setdiff(trainingGroups, eg);

        trIndex = findIdx(tg, stimRespData.groupIndex);
        esIndex = findIdx(eg, stimRespData.groupIndex);    
        tic;
        [modelParamsES, optOptionsES] = strfOpt(modelParams, trIndex, optOptions, esIndex);
        etime = toc;
        fprintf('\nEarly stopping Run took %0.0fm\n', etime/60);
 
        strfs(k, :, :) = squeeze(modelParamsES.w1);
        biases(k) = modelParamsES.b1;
        numIters(k) = optOptionsES.diagnostics.bestiter;
        runTimes(k) = etime;
        if strcmp(modelType, 'leglm')
            ms(k) = modelParamsES.m;
        end
        
        totalTime = totalTime + etime;
        
        clear modelParamsES;
        clear optOptionsES;
    end
    
    %% compute the final STRF and bias
    finalStrf = squeeze(mean(strfs, 1));
    size(finalStrf)
    finalBias = mean(biases);
    if strcmp(modelType, 'leglm')
        finalM = mean(ms);
    end
    
    modelParamsTrained = modelParams;
    modelParamsTrained.w1 = finalStrf;
    modelParamsTrained.b1 = finalBias;
    if strcmp(modelType, 'leglm')
        modelParamsTrained.m = finalM;
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
    
    trainingIndex = findIdx(trainingGroups, stimRespData.groupIndex);
    validationIndex = findIdx(validationGroups, stimRespData.groupIndex);
    
    %% produce responses
    [modelResp, rawModelResp] = compute_response(stimRespData, modelParamsTrained);
    if strcmp(modelType, 'poisson')
        modelResp = rawModelResp;
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
    h5.set_ds(fid, '/data', 'training_groups', trainingGroups);
    h5.set_ds(fid, '/data', 'validation_groups', validationGroups);
    
    h5.set_ds(fid, '/data', 'holdout_subgroups', earlyStoppingGroups);
    
    h5.set_attr(fid, '/opt', 'method', 'threshgrad');
    h5.set_ds(fid, '/opt', 'threshold', optOptions.threshold);
    h5.set_ds(fid, '/opt', 'max_iter', optOptions.maxIter);
    h5.set_ds(fid, '/opt', 'num_iters', numIters);
    h5.set_ds(fid, '/opt', 'run_times', runTimes);
    h5.set_ds(fid, '/opt', 'final_run_time', totalTime);
    h5.set_ds(fid, '/opt', 'version', 5);
            
    h5.set_attr(fid, '/model', 'strf_length', strfLength);
    h5.set_attr(fid, '/model', 'transforms', trStrs);
    h5.set_attr(fid, '/model', 'type', modelType);
    h5.set_ds(fid, '/model', 'cv_weights', strfs);
    h5.set_ds(fid, '/model', 'cv_bias', biases);
    h5.set_ds(fid, '/model', 'weights', finalStrf);
    h5.set_ds(fid, '/model', 'bias', finalBias);
    h5.set_ds(fid, '/model', 'response', modelResp);    
    h5.set_ds(fid, '/model', 'response_raw', rawModelResp);    
    if strcmp(modelType, 'leglm')
        h5.set_ds(fid, '/model', 'm', finalM);
    end
    
    h5.close(fid);
    
end

function modelParams = init_glm_model(stimRespData, strfDelays, modelType)

    switch modelType
        case 'gaussian'
            family = init_glm_family_gaussian();
            outputNL = @(x) identity_outputnl(x);            
        case 'binomial'
            numTrials = round(mean(stimRespData.numTrials));
            family = init_glm_family_binomial(numTrials);
            outputNL = @(x) sigmoid_outputnl(x);
        case 'poisson'
            family = init_glm_family_poisson();
            outputNL = @(x) exp_outputnl(x);
        case 'leglm'
            family = init_glm_family_poisson();            
    end
    
    dispersion = 1;
    if ~strcmp(modelType, 'leglm')
        modelParams = glmInit(stimRespData.numChannels, strfDelays, family, outputNL, dispersion);
    else
        modelParams = leglmInit(stimRespData.numChannels, strfDelays, family, dispersion);
    end
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
