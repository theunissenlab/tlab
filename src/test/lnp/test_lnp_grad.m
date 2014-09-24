function test_lnp_grad

    nDelays = 3;
    nChannels = 2;
    delays = 0:(nDelays-1);
    modelParams = lnpInit(nChannels, delays);
    
    %w = [0.25 0.4];
    w = [0.25 0.4 0.3; 0.1 0.2 0.5];
    B = 0.25;
    
    modelParams.w1 = w;
    modelParams.b1 = 0.25;
    modelParams.regularize = 1;
    modelParams.regularizeWeight = 1;
    
    spikeTimes = zeros(100, 1);
    spikeIndicies = [5 24 36 55 80 90];
    spikeTimes(spikeIndicies) = 1;
    
    stim = randn(100, nChannels);
    stim(spikeIndicies) = 30;
    
    spikeTimes(24) = 2;
    spikeTimes(55) = 4;
    
    groupIndex = ones(100, 1);
    
    global globDat;
    strfData(stim, spikeTimes, groupIndex);
    
    [modelParams, gfd] = lnpGradFD(modelParams, []);
    gfd
    
    [modelParams, g] = lnpGrad(modelParams, []);
    g
    
    