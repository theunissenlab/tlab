function test_lnpsr2_simple2()

    sampleRate = 1000;
    duration = 2;
    nDataPoints = round(duration * sampleRate);

    nStimDelays = 25;
    numStimChannels = 10;
    stimDelays = 0:(nStimDelays-1);
    spikeDelays = [];
    
    nlType = 'exponential';
        
    %% initialize model
    ignoreBias = 1;
    modelParamsReal = lnpsr2Init(numStimChannels, stimDelays, nlType, ...
				 spikeDelays, sampleRate, ignoreBias);

    %w = [0.75 0.4]*3;
    w = randn(numStimChannels, nStimDelays);
    w = (w / norm(w(:))) * 3;
    
    modelParamsReal.w1 = w;
    modelParamsReal.b1 = 0;
    modelParamsReal.spikeResponseWeights = [];
    
    %% create test data
    stim = randn(nDataPoints, numStimChannels);
    
    % set up global structure
    groupIndex = ones(nDataPoints, 1);
    global globDat;
    strfData(stim, zeros(1, nDataPoints), groupIndex);    

    numRealTrials = 20;    
    datIdx = 1:nDataPoints;
    [modelParamsReal, modelResponseReal, fullResponse] = lnpsr2Fwd(modelParamsReal, datIdx, numRealTrials);
    
    realSpikeTrials = fullResponse.spikeTrials(1, :);
    realSpikeTimes = (find(realSpikeTrials > 0) - 1) / sampleRate;
    sint = 1/sampleRate;
    t = 0:sint:(nDataPoints-1)*sint;
    
    figure; hold on;
    subplot(4, 1, 1);
    plot(t, stim, 'k-');
    title('Stimulus');
    subplot(4, 1, 2);
    plot(t, fullResponse.stimCurrent, 'b-');
    title('Linear Response');
    subplot(4, 1, 3); hold on;
    exampleSets = [1:10:numRealTrials];
    for k = 1:length(exampleSets)
      plot(t, fullResponse.nonlinearResponse(exampleSets(k), :));
    end
    title('Example Poisson Rates');
    subplot(4, 1, 4); hold on;
    plot(t, modelResponseReal, 'k-');
    title('PSTH');
        
    %% Set options for optimization method
    optOptions = trnSCG;
    %optOptions.minErrChange = 1e-2;
    %optOptions.minStepSize = 1e-2;
    optOptions.display = 1;
        
    %% set up global structure with fake data
    groupIndex = ones(nDataPoints, 1);
    
    spikeTimesCell = cell(numRealTrials, 1);
    for k = 1:numRealTrials
      strials = fullResponse.spikeTrials(k, :);
      stimes = (find(strials > 0) - 1) / sampleRate;
      spikeTimesCell{k} = stimes;      
    end
    strfData(stim, spikeTimesCell, groupIndex);
    
    %% Initialize LNP spike-response model
    modelParams = lnpsr2Init(numStimChannels, stimDelays, nlType, ...
			     spikeDelays, sampleRate, modelParamsReal.ignoreBias);
    
    fprintf('Using random initial guess...\n');
    modelParams.w1 = randn(size(modelParamsReal.w1));
    modelParams.b1 = 0;
    
    modelParams.checkGrad = 0;    
    modelParams.useFDGrad = 0;
    
    %% run strflab optimization
    [modelParamsTrained, optOptions] = strfOpt(modelParams, datIdx, optOptions);
       
    actualParams = modelParamsReal.w1
    fitParams = modelParamsTrained.w1
 
    %% look at response
    numTrials = 100;
    [modelParamsTrained, predResp] = lnpsr2Fwd(modelParamsTrained, datIdx, numTrials);
    predResp = (predResp/max(predResp)) * max(modelResponseReal);
    
    %% compute coherence    
    spikeTimesMs = cell(length(spikeTimesCell), 1);
    for k = 1:length(spikeTimesCell)
        spikeTimesMs{k} = spikeTimesCell{k} * 1e3;
    end
    stimLengthMs = (nDataPoints / sampleRate) * 1e3;    
    psthdata = split_psth(spikeTimesMs, stimLengthMs);
    psthHalf1 = psthdata.psth_half1;
    psthHalf2 = psthdata.psth_half2;
    
    infoFreqCutoff = 90;
    infoWindowSize = 0.500;
    [cBound, cModel] = compute_coherence_full(predResp, modelResponseReal, psthHalf1, psthHalf2, sampleRate, numRealTrials, infoFreqCutoff, infoWindowSize);
        
    performanceRatio = cModel.info / cBound.info
        
    realStrf = modelParamsReal.w1
    predStrf = modelParamsTrained.w1
    strfDiff = norm(realStrf - predStrf)
  
    realw0 = modelParamsReal.b1
    predw0 = modelParamsTrained.b1
    
    realSwts = modelParamsReal.spikeResponseWeights
    predSwts = modelParamsTrained.spikeResponseWeights
    
    figure; hold on;
    plot(modelResponseReal, 'k-', 'LineWidth', 2);
    plot(predResp, 'r-', 'LineWidth', 2);
    axis tight;
    %axis([1 length(modelResponseReal) 0 1]);

    %% plot information
    figure; hold on;
    plot(cBound.f, cBound.c, 'k-', 'LineWidth', 2);
    plot(cBound.f, cBound.cUpper, 'b-', 'LineWidth', 2);
    plot(cBound.f, cBound.cLower, 'r-', 'LineWidth', 2);

    plot(cModel.f, cModel.c, 'k--', 'LineWidth', 2);
    plot(cModel.f, cModel.cUpper, 'b--', 'LineWidth', 2);
    plot(cModel.f, cModel.cLower, 'r--', 'LineWidth', 2);
    theTitle = sprintf('Info=%0.0f out of %0.0f bits | Ratio=%0.2f', ...
                       cModel.info, cBound.info, performanceRatio);
    title(theTitle);
    axis([min(cBound.f), max(cBound.f), 0, 1]);
    
    %save('~/test_lnpsr2_data.mat', 'stim', 'predResp', 'modelResponseReal', 'modelParamsTrained', 'modelParamsReal', 'cBound', 'cModel', 'spikeTimes');
%}