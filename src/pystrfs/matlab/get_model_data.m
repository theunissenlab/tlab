function mdata = get_model_data(responseFile)

    h5 = h5utils();
    fid = h5.open(responseFile);
    
    unitFile = h5.get_attr(fid, '/data', 'unit_file');
    preprocFile = h5.get_attr(fid, '/data', 'preproc_file');
    
    trStrs = h5.get_attr(fid, '/model', 'transforms');
    trStrs = strtrim(trStrs);
    if ~isempty(trStrs)
        transforms = regexp(trStrs,',','split');
    else
        transforms = '';
    end
    
    data = struct;
    data.stimClass = h5.get_attr(fid, '/data', 'stim_class');
    data.sampleRate = h5.get_attr(fid, '/data', 'sample_rate');
    data.groupIndex = h5.get_ds(fid, '/data/group_index');
    data.trainingIndex = h5.get_ds(fid, '/data/training_index');
    data.validationIndex = h5.get_ds(fid, '/data/validation_index');
    
    stimRespData = get_stimresp_data(unitFile, preprocFile, data.stimClass);    
    data.f = stimRespData.f;
    data.numChannels = stimRespData.numChannels;    
    data.allSpikeTrials = stimRespData.allSpikeTrials;
    data.wholeResp = stimRespData.wholeResp;
    data.numTrials = stimRespData.numTrials;
    
    wholeStim = stimRespData.wholeStim;
    for k = 1:length(transforms)           
        tfunc = sprintf('transform_%s(wholeStim)', transforms{k});
        wholeStim = eval(tfunc);        
    end    
    data.wholeStim = wholeStim;    
    
    model = struct;
    model.response = h5.get_ds(fid, '/model/response');    
    model.rawResponse = h5.get_ds(fid, '/model/response_raw');
    model.transforms = transforms;
    
    modelType = h5.get_attr(fid, '/model', 'type');
    optMethod = h5.get_attr(fid, '/opt', 'method');
    
    model.type = modelType;
    
    opt = struct;
    opt.optMethod = optMethod;
    
    %% get STRF    
    if ~ismember(modelType, {'nl_dists', 'nl_spline'})        
        model.strf = h5.get_ds(fid, '/model/weights');
        model.bias = h5.get_ds(fid, '/model/bias');
    else
        origFile = h5.get_attr(fid, '/model', 'original_response_file');
        fid2 = h5.open(origFile);
        model.strf = h5.get_ds(fid2, '/model/weights');
        model.bias = h5.get_ds(fid2, '/model/bias');
        h5.close(fid2);
    end
        
    %% get output NL
    if ismember(modelType, {'nl_dists', 'nl_spline'})
        [rootDir, baseName, ext] = fileparts(responseFile);
        matFile = fullfile(rootDir, [baseName '.mat']);
        vars = load(matFile);
        model.nlinfo = vars.nlinfo;
        model.outputNL = @(x) output_nl_wrapper(x, vars.nlinfo.outputNL);
        clear vars;
    else       
        switch model.type            
            case 'poisson'
                model.outputNL = @(x) exp(x);                
            case 'binomial'
                model.outputNL = @(x) 1 ./ (1 + exp(-x));                
            case 'leglm'
                model.m = h5.get_ds(fid, '/model/m');
                model.outputNL = @(x) log(1 + exp(x)).^model.m;                
            case 'linear'
                model.outputNL = @(x) x;
        end
    end
    
    
    coherence = struct;
    coherence.bound.training.upper = h5.get_ds(fid, '/model/performance/coherence/bound/training/coherence_upper');
    coherence.bound.training.lower = h5.get_ds(fid, '/model/performance/coherence/bound/training/coherence_lower');
    coherence.bound.training.mean = h5.get_ds(fid, '/model/performance/coherence/bound/training/coherence_mean');
    coherence.bound.training.f = h5.get_ds(fid, '/model/performance/coherence/bound/training/frequency');
    coherence.bound.training.info_upper = h5.get_ds(fid, '/model/performance/coherence/bound/training/info_upper');
    coherence.bound.training.info_lower = h5.get_ds(fid, '/model/performance/coherence/bound/training/info_lower');
    coherence.bound.training.info_mean = h5.get_ds(fid, '/model/performance/coherence/bound/training/info_mean');
    
    coherence.bound.validation.upper = h5.get_ds(fid, '/model/performance/coherence/bound/validation/coherence_upper');
    coherence.bound.validation.lower = h5.get_ds(fid, '/model/performance/coherence/bound/validation/coherence_lower');
    coherence.bound.validation.mean = h5.get_ds(fid, '/model/performance/coherence/bound/validation/coherence_mean');
    coherence.bound.validation.f = h5.get_ds(fid, '/model/performance/coherence/bound/validation/frequency');
    coherence.bound.validation.info_upper = h5.get_ds(fid, '/model/performance/coherence/bound/validation/info_upper');
    coherence.bound.validation.info_lower = h5.get_ds(fid, '/model/performance/coherence/bound/validation/info_lower');
    coherence.bound.validation.info_mean = h5.get_ds(fid, '/model/performance/coherence/bound/validation/info_mean');
    
    coherence.training.upper = h5.get_ds(fid, '/model/performance/coherence/training/coherence_upper');
    coherence.training.lower = h5.get_ds(fid, '/model/performance/coherence/training/coherence_lower');
    coherence.training.mean = h5.get_ds(fid, '/model/performance/coherence/training/coherence_mean');
    coherence.training.f = h5.get_ds(fid, '/model/performance/coherence/training/frequency');
    coherence.training.info_upper = h5.get_ds(fid, '/model/performance/coherence/training/info_upper');
    coherence.training.info_lower = h5.get_ds(fid, '/model/performance/coherence/training/info_lower');
    coherence.training.info_mean = h5.get_ds(fid, '/model/performance/coherence/training/info_mean');
    
    coherence.validation.upper = h5.get_ds(fid, '/model/performance/coherence/validation/coherence_upper');
    coherence.validation.lower = h5.get_ds(fid, '/model/performance/coherence/validation/coherence_lower');
    coherence.validation.mean = h5.get_ds(fid, '/model/performance/coherence/validation/coherence_mean');
    coherence.validation.f = h5.get_ds(fid, '/model/performance/coherence/validation/frequency');
    coherence.validation.info_upper = h5.get_ds(fid, '/model/performance/coherence/validation/info_upper');
    coherence.validation.info_lower = h5.get_ds(fid, '/model/performance/coherence/validation/info_lower');
    coherence.validation.info_mean = h5.get_ds(fid, '/model/performance/coherence/validation/info_mean');
    
    h5.close(fid);    
    
    [rootDir, baseName, ext] = fileparts(responseFile);
    spath = regexp(rootDir, '/', 'split');    
    cellName = spath{end-1};
    
    mdata = struct;
    mdata.cellName = cellName;
    mdata.name = baseName;
    mdata.responseFile = responseFile;
    mdata.data = data;
    mdata.model = model;
    mdata.coherence = coherence;
    mdata.opt = opt;
    
end

function y = output_nl_wrapper(x, fnObj)    
    y = fnval(fnObj, x);
end
