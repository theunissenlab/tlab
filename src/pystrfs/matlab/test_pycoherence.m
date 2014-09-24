function test_pycoherence(cfile)

    %% get data
    h5 = h5utils();
    f = h5.open(cfile, 'r');
    
    tnums = h5.get_subgroups(f, '/trials');
    numTrials = length(tnums);
    
    allTrials = cell(numTrials, 1);
    for k = 1:numTrials
        ds = sprintf('/trials/%s', tnums{k});
        trial = h5.get_ds(f, ds);
        allTrials{k} = trial;        
    end
    allTrials = {allTrials};
    
    psthModel = h5.get_ds(f, '/psth_model');
    psthModel = psthModel';
    
    cpyBoundF = h5.get_ds(f, '/cpy/bound/frequency');
    cpyBoundMean = h5.get_ds(f, '/cpy/bound/coherence');
    
    cpyModelF = h5.get_ds(f, '/cpy/model/frequency');
    cpyModelMean = h5.get_ds(f, '/cpy/model/coherence');
    
    cpyNoNormF = h5.get_ds(f, '/cpy/nonorm/frequency');
    cpyNoNormMean = h5.get_ds(f, '/cpy/nonorm/coherence');
    
    binSize = 0.001;
    sampleRate = 1 / binSize;
    duration = length(psthModel) * binSize;
            
    %% compute PSTH and split PSTHs
    preprocRespParams = struct;         %create preprocessing struct
    preprocRespParams.units = 's';      %files have units of s
    preprocRespParams.binSize = binSize;  %1ms bin size, same sample rate as spectrograms
    preprocRespParams.stimLengths = [duration]; %needed to compute correct PSTH lengths

    [wholeResp, groupIndex, respInfo, preprocRespParams] = preprocSpikeResponses(allTrials, preprocRespParams);

    preprocRespParams.split = 1;
    [wholeSplitResp, respSplitInfo, preprocRespParams] = preprocSpikeResponses(allTrials, preprocRespParams);

    infoFreqCutoff = -1;
    infoWindowSize = 0.200;
        
    %% Compute model responses, coherences, and un-zscore them
    [cBound, cModel] = compute_coherence_full(psthModel, wholeResp, wholeSplitResp(1, :), wholeSplitResp(2, :), sampleRate, numTrials, infoFreqCutoff, infoWindowSize);
    
    %% compute non-normed model
    cNoNorm = compute_coherence_mean(psthModel, wholeResp, sampleRate, infoFreqCutoff, infoWindowSize);
    
    %% plot
    figure(); hold on;    
    subplot(4, 1, 1);
    plot(wholeResp, 'k-');
    subplot(4, 1, 2);
    plot(psthModel, 'r-');
    subplot(4, 1, 3);
    plot(psthModel - wholeResp, 'g-');
    
    subplot(4, 1, 4); hold on;
    plot(cBound.f, cBound.c, 'k-');
    plot(cModel.f, cModel.c, 'r-');
    plot(cpyBoundF, cpyBoundMean, 'k--');
    plot(cpyModelF, cpyModelMean, 'r--');
    
    title(sprintf('Info: Bound=%.2f (%d Hz), Model=%.2f (%d Hz), Ratio=%.3f', cBound.info, cBound.freqCutoff, cModel.info, cModel.freqCutoff, cModel.info / cBound.info));
    
    figure(); hold on;
    plot(cNoNorm.f, cNoNorm.c, 'k-');
    plot(cpyNoNormF, cpyNoNormMean, 'r-');
    legend('MATLAB', 'Python');
    
    h5.close(f);

