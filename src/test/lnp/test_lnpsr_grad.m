function test_lnpsr_grad

    nDataPoints = 100;

    nStimDelays = 3;
    numStimChannels = 2;
    stimDelays = 0:(nStimDelays-1);
        
    nSpikeDelays = 4;
    spikeDelays = 0:(nSpikeDelays-1);
    
    nlType = 'exponential';
    %nlType = 'logexp1';
    
    sampleRate = 1000;

    %% initialize model
    numTrials = 50;
    modelParamsReal = lnpsrInit(numStimChannels, stimDelays, nlType, spikeDelays, sampleRate, numTrials);
    
    %w = [0.25 0.4 0.3; 0.1 0.2 0.5];
    w = [0.25 0.10 0.05; 0.15 0.05 0.01];
    
    modelParamsReal.w1 = w;
    modelParamsReal.b1 = 0.0;
    modelParamsReal.spikeResponseWeights = [-1.5 -1.5 -1 -0.5];
    
    %% create test data
    gain = 4;
    stim = gen_test_stim(nDataPoints, numStimChannels, sampleRate)*gain;

    % set up global structure
    groupIndex = ones(nDataPoints, 1);
    global globDat;
    strfData(stim, zeros(1, nDataPoints), groupIndex);

    datIdx = 1:nDataPoints;
    [modelParamsReal, modelResponseReal] = lnpsrFwd(modelParamsReal, datIdx, numTrials);
    strfData(stim, modelResponseReal, groupIndex);
    
    modelParamsReal.l2weight = 1;
    
    modelParamsReal.regularize = 0;
    modelParamsReal.regularizeWeight = 1e6;
    
    [modelParamsReal, gfd] = lnpsrGradFD(modelParamsReal, datIdx);
    gfd
    
    [modelParamsReal, g] = lnpsrGrad(modelParamsReal, datIdx);
    g
    