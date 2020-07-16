function test_lnpsr2

    %nDataPoints = 1000;
    sampleRate = 1000;
    nDataPoints = 5000;

    %nStimDelays = 3;
    %numStimChannels = 2;
    nStimDelays = 2;
    numStimChannels = 1;
    stimDelays = 0:(nStimDelays-1);
    
    useSrFilt = 0;
    if useSrFilt
        nSpikeDelays = 4;
        spikeDelays = 0:(nSpikeDelays-1);
    else
        nSpikeDelays = 0;
        spikeDelays = [];
    end
    
    
    nlType = 'exponential';
    %nlType = 'logexp1';
        
    %% initialize model
    
    modelParamsReal = lnpsr2Init(numStimChannels, stimDelays, nlType, spikeDelays, sampleRate);
    
    %w = [1 0.65 0.05; 0.8 0.7 0.1];
    w = [0.75 0.4]*2;
    b = 0.0;    
    
    modelParamsReal.w1 = w;
    modelParamsReal.b1 = b;
    modelParamsReal.ignoreBias = 1;
    if useSrFilt
        %modelParamsReal.spikeResponseWeights = [-1.5 -1.5 -1 -0.5];    
        %modelParamsReal.spikeResponseWeights = [1.5 1 0.5 0.25 0 -0.5 -0.5 -0.25 -0.1];    
        modelParamsReal.spikeResponseWeights = [-1 -0.75 -0.50 -0.25];
        %modelParamsReal.spikeResponseWeights = [-2];
    else
        modelParamsReal.spikeResponseWeights = [];
    end
    
    %% create test data
    gain = 2;
    stim = randn(nDataPoints, numStimChannels);
    %stim = gen_test_stim(nDataPoints, numStimChannels, sampleRate)*gain;
    
    %subtract off mean stimulus
    %stim = stim - mean(stim);
    
    % set up global structure
    groupIndex = ones(nDataPoints, 1);
    global globDat;
    strfData(stim, zeros(1, nDataPoints), groupIndex);
    
    numTrials = 100;
    numRealTrials = 1000;
    
    datIdx = 1:nDataPoints;
    [modelParamsReal, modelResponseReal, fullResponse] = lnpsr2Fwd(modelParamsReal, datIdx, numRealTrials);

    %subtract of mean response
    %modelResponseReal = modelResponseReal - mean(modelResponseReal);
    
    sint = 1/sampleRate;
    t = 0:sint:(nDataPoints-1)*sint;
    
    makePlots = 0;
    if makePlots
    
        for k = 1:numTrials
        
            h1 = figure; hold on;
           
            %subplot(2, 1, 1); hold on;
            plot(t, fullResponse.stimCurrent, 'b-');
            plot(t, fullResponse.spikeCurrent(k, :), 'r-');
            legend('Stim', 'Spike');
            title(sprintf('%d: Linear Response', k));
            axis tight;
            plot2svg_2d('~/temp/lnpsr2-linresp.svg', h1);

            h2 = figure; hold on;
            spIndx = find(fullResponse.spikeTrials(k, :) > 0);
            %subplot(2, 1, 2); hold on;
            plot(t, fullResponse.nonlinearResponse(k, :), 'b-');
            plot(t(spIndx), ones(1, length(spIndx)), 'ro');
            title(sprintf('%d: Nonlinear Response', k));
            %axis([min(t) max(t) 0 5]);
            axis tight;
            plot2svg_2d('~/temp/lnpsr2-nonlinresp.svg', h2);
        end
    end
    
    figure; hold on;
    subplot(2, 1, 1);
    imagesc(t, size(stim, 2), stim');
    colorbar;
    title('Stimulus');

    subplot(2, 1, 2);
    plot(t, modelResponseReal, 'k-');
    title('Model Response');

    
    %% Set options for optimization method
    %{
    optOptions = trnGradDescLineSearch;
    optOptions.display = 1;
    optOptions.coorDesc = 0;
    optOptions.earlyStop = 0;
    optOptions.stepSize = 1e-4;
    optOptions.nDispTruncate = 0;
    optOptions.maxIter = 15000;
    %optOptions.lineSearch = 1;
    %}
    
    %optOptions = trnSCGNetlab;
    optOptions = trnSCG;
    %optOptions.minErrChange = 1e-2;
    %optOptions.minStepSize = 1e-2;
    
    optOptions.display = 1;
        
    %% set up global structure with fake data
    groupIndex = ones(nDataPoints, 1);
    
    spikeTimes = cell(numRealTrials, 1);
    for k = 1:numRealTrials
        st = fullResponse.spikeTrials(k, :);
        spikeTimes{k} = (find(st > 0) - 1) / sampleRate;
    end
    strfData(stim, spikeTimes, groupIndex);
    
    %% Initialize LNP spike-response model
    nlType = 'exponential';
    modelParams = lnpsr2Init(numStimChannels, stimDelays, nlType, ...
			     spikeDelays, sampleRate);
    modelParams.ignoreBias = 1;
    
    modelParams.fixedStrf = 0;
    if modelParams.fixedStrf
        fprintf('Using fixed STRF...\n');
        modelParams.w1 = modelParamsReal.w1;
        modelParams.b1 = modelParamsReal.b1;
    end 
        
    usePerturbedInitialGuess = 0;
    if usePerturbedInitialGuess
        fprintf('Using perturbed initial guess...\n');
        rgain = 1e-3;
        strf = modelParamsReal.w1;
        modelParams.w1 = strf + randn(size(strf))*mean(mean(strf))*rgain;
        modelParams.b1 = modelParamsReal.b1 + randn()*modelParamsReal.b1*rgain;
        swts = modelParams.spikeResponseWeights;
        modelParams.spikeResponseWeights = swts + randn(size(swts))*mean(swts)*rgain;
    end
    
    useRandomInitialGuess = 1;
    if useRandomInitialGuess
        fprintf('Using random initial guess...\n');
        rgain = 1e-6;
        strf = modelParamsReal.w1;
        modelParams.w1 = randn(size(strf));
        modelParams.b1 = 0;
        swts = modelParamsReal.spikeResponseWeights;
        modelParams.spikeResponseWeights = randn(size(swts))*mean(swts)*rgain;
    end
    
    modelParams.checkGrad = 0;    
    modelParams.useFDGrad = 0;
    
    %% run strflab optimization
    [modelParamsTrained, optOptions] = strfOpt(modelParams, datIdx, optOptions);
       
    %% look at response
    [modelParamsTrained, predResp] = lnpsr2Fwd(modelParamsTrained, datIdx, numTrials);
    predResp = (predResp/max(predResp)) * max(modelResponseReal);
    
    %% compute coherence    
    spikeTimesMs = cell(length(spikeTimes), 1);
    for k = 1:length(spikeTimes)
        spikeTimesMs{k} = spikeTimes{k} * 1e3;
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
    
    save('~/test_lnpsr2_data.mat', 'stim', 'predResp', 'modelResponseReal', 'modelParamsTrained', 'modelParamsReal', 'cBound', 'cModel', 'spikeTimes');
    