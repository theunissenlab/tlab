function test_lnpsr2_simple()

    sampleRate = 1000;
    duration = 500;
    nDataPoints = round(duration * sampleRate);

    nStimDelays = 1;
    numStimChannels = 1;
    stimDelays = 0:(nStimDelays-1);

    nSpikeDelays = 3;
    spikeDelays = 0:(nSpikeDelays-1);
    
    nlType = 'exponential';
        
    %% initialize model
    ignoreBias = 1;
    modelParamsReal = lnpsr2Init(numStimChannels, stimDelays, nlType, ...
				 spikeDelays, sampleRate, ignoreBias);

    %w = [0.75 0.4]*2;
    w = [1];
    %s = [4 2.5 1.5 0.5];
    s = [3 2 1];
    %w = [0.75 0.4 -0.2; 0.6 0.45 -0.1];
    
    modelParamsReal.w1 = w;
    modelParamsReal.b1 = 0;
    modelParamsReal.spikeResponseWeights = s;
    
    %% create test data
    stim = randn(nDataPoints, numStimChannels) + 1;
    
    % set up global structure
    groupIndex = ones(nDataPoints, 1);
    global globDat;
    strfData(stim, zeros(1, nDataPoints), groupIndex);    

    numRealTrials = 1;    
    datIdx = 1:nDataPoints;
    [modelParamsReal, modelResponseReal, fullResponse] = lnpsr2Fwd(modelParamsReal, datIdx, numRealTrials);
    
    realSpikeTrials = fullResponse.spikeTrials(1, :);
    realSpikeTimes = (find(realSpikeTrials > 0) - 1) / sampleRate;
    fprintf('# of event times: %d\n', length(realSpikeTimes));
    sint = 1/sampleRate;
    t = 0:sint:(nDataPoints-1)*sint;
    
    figure; hold on;
    subplot(3, 1, 1);
    plot(t, stim, 'k-');
    title('Stimulus');
    subplot(3, 1, 2);
    plot(t, fullResponse.stimCurrent, 'b-');
    title('Linear Response');
    subplot(3, 1, 3); hold on;
    plot(t, fullResponse.nonlinearResponse(1, :), 'g-');
    plot(realSpikeTimes, zeros(size(realSpikeTimes)), 'r.');
    title('Poisson Rate');
        
    %% Set options for optimization method
    optOptions = trnSCG;
    %optOptions.minErrChange = 1e-2;
    %optOptions.minStepSize = 1e-2;
    optOptions.display = 1;
        
    %% set up global structure with fake data
    groupIndex = ones(nDataPoints, 1);
    
    spikeTimesCell = {realSpikeTimes}; 
    strfData(stim, spikeTimesCell, groupIndex);
    
    %% Initialize LNP spike-response model
    modelParams = lnpsr2Init(numStimChannels, stimDelays, nlType, ...
			     spikeDelays, sampleRate, modelParamsReal.ignoreBias);
    
    %{
    fprintf('Using random initial guess...\n');
    modelParams.w1 = randn(size(modelParamsReal.w1));
    modelParams.b1 = 0;
    %}
    modelParams.fixedStrf = 1;
    modelParams.w1 = modelParamsReal.w1;
    modelParams.spikeResponseWeights = randn(size(modelParamsReal.spikeResponseWeights)) - 1;
    siguess = modelParams.spikeResponseWeights
    
    modelParams.checkGrad = 0;    
    modelParams.useFDGrad = 0;
    modelParams.useLowMem = 1;
    
    %% run strflab optimization
    [modelParamsTrained, optOptions] = strfOpt(modelParams, datIdx, optOptions);
       
    actualStrf = modelParamsReal.w1
    fitStrf = modelParamsTrained.w1
    
    actualSrWts = modelParamsReal.spikeResponseWeights
    fitSrWts = modelParamsTrained.spikeResponseWeights
    
    