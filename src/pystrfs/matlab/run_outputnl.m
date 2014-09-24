function run_outputnl(responseFile, outputFile, nlType)
    
    %% get parameters from optimization output file
    h5 = h5utils();
    fid = h5.open(responseFile);
    
    unitFile = h5.get_attr(fid, '/data', 'unit_file');
    unit = read_unit_file(unitFile);
    
    preprocFile = h5.get_attr(fid, '/data', 'preproc_file');
    
    stimClass = h5.get_attr(fid, '/data', 'stim_class');
    sampleRate = h5.get_attr(fid, '/data', 'sample_rate');
    
    groupIndex = h5.get_ds(fid, '/data/group_index');
    trainingIndex = h5.get_ds(fid, '/data/training_index');
    validationIndex = h5.get_ds(fid, '/data/validation_index');
    trainingGroups = h5.get_ds(fid, '/data/training_groups');
    validationGroups = h5.get_ds(fid, '/data/validation_groups');
        
    rawModelResp = h5.get_ds(fid, '/model/response_raw');
    trStrs = h5.get_attr(fid, '/model', 'transforms');
    
    optMethod = h5.get_attr(fid, '/opt', 'method');
    
    
    %% get spike data from unit
    allResps = unit.class_responses.(stimClass);
    allSpikeTrials = cell(length(allResps), 1);
    numTrials = zeros(length(allResps), 1);
    respLengths = zeros(length(allResps), 1);
    for k = 1:length(allResps)       
        allSpikeTrials{k} = allResps{k}.spikeTrials;        
        numTrials(k) = allResps{k}.numTrials;
        respLengths(k) = sum(groupIndex == k) / sampleRate;
    end
    
    
    %% compute PSTH and split PSTHs
    preprocRespParams = struct;         %create preprocessing struct
    preprocRespParams.units = 's';      %files have units of s
    preprocRespParams.binSize = 0.001;  %1ms bin size, same sample rate as spectrograms
    preprocRespParams.stimLengths = respLengths; %needed to compute correct PSTH lengths

    [wholeResp, groupIndex, respInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);
    
    avgNumTrials = mean(numTrials);
    
    if strcmp('dists', nlType)
        numResample = 1000;
        resampleFraction = 0.75;
        nlinfo = estimate_spike_dists_cv(cv(rawModelResp(trainingIndex)), cv(wholeResp(trainingIndex)), avgNumTrials, numResample, resampleFraction);
        %nlinfo = estimate_spike_dists_cv2(cv(rawModelResp(trainingIndex)), cv(wholeResp(trainingIndex)), avgNumTrials, numDensityPoints, numResample, resampleFraction);
    elseif strcmp('spline', nlType)
        numResample = 75;
        resampleFraction = 0.75;
        nlinfo = estimate_outputnl_from_spline_cv(cv(rawModelResp(trainingIndex)), cv(wholeResp(trainingIndex)), numResample, resampleFraction);
    end
    
    nlModelResp = fnval(nlinfo.outputNL, rawModelResp);
    tname = sprintf('nl_%s', nlType);   
    
    [rootDir, baseName, ext] = fileparts(outputFile);
    outputFileMat = [rootDir '/' baseName '.mat'];
    
    if exist(outputFile, 'file')
        delete(outputFile);
    end
    if exist(outputFileMat, 'file')
        delete(outputFileMat);
    end
    
    save(outputFileMat, 'nlinfo');
    
    h5 = h5utils();
    fid = h5.create(outputFile);    
    
    h5.set_attr(fid, '/data', 'sample_rate', sampleRate);
    h5.set_attr(fid, '/data', 'unit_file', unitFile);
    h5.set_attr(fid, '/data', 'preproc_file', preprocFile);
    h5.set_attr(fid, '/data', 'stim_class', stimClass);    
    
    h5.set_ds(fid, '/data', 'training_groups', trainingGroups);
    h5.set_ds(fid, '/data', 'validation_groups', validationGroups);
    h5.set_ds(fid, '/data', 'group_index', groupIndex);
    h5.set_ds(fid, '/data', 'training_index', trainingIndex);
    h5.set_ds(fid, '/data', 'validation_index', validationIndex);
    
    h5.set_attr(fid, '/model', 'type', tname);
    h5.set_attr(fid, '/model', 'transforms', trStrs);
    h5.set_attr(fid, '/model', 'original_response_file', responseFile);
    h5.set_attr(fid, '/model', 'mat_file', outputFileMat);
    
    h5.set_attr(fid, '/opt', 'method', optMethod);    
    
    h5.set_ds(fid, '/model', 'response', nlModelResp);
    h5.set_ds(fid, '/model', 'response_raw', nlModelResp);
    
    h5.close(fid);
    