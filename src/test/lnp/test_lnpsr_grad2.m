function test_lnpsr_grad2

    duration = 1;
    sampleRate = 1000;
    nDataPoints = round(duration*sampleRate);

    nStimDelays = 2;
    numStimChannels = 2;
    stimDelays = 0:(nStimDelays-1);
        
    nSpikeDelays = 2;
    spikeDelays = 0:(nSpikeDelays-1);
    
    nlType = 'exponential';
    %nlType = 'harris';

    %% initialize model
    ignoreBias = 0;
    numTrials = 50;    
    modelParamsReal = lnpsr2Init(numStimChannels, stimDelays, nlType, spikeDelays, sampleRate, ignoreBias);
    
    modelParamsReal.w1 = randn(numStimChannels, nStimDelays);
    modelParamsReal.b1 = 0.5;
    modelParamsReal.spikeResponseWeights = randn(1, nSpikeDelays);
    
    %% create test data
    stim = randn(nDataPoints, numStimChannels);

    % set up global structure
    groupIndex = ones(nDataPoints, 1);
    global globDat;
    strfData(stim, zeros(1, nDataPoints), groupIndex);

    datIdx = 1:nDataPoints;
    [modelParamsReal, modelResponseReal, fullResponse] = lnpsr2Fwd(modelParamsReal, datIdx, numTrials);
    
    %convert response spike times to seconds
    spikeTimes = cell(numTrials, 1);
    for k = 1:numTrials
        st = fullResponse.spikeTrials(k, :);
        spikeTimes{k} = (find(st > 0) - 1) / sampleRate;
    end
    
    strfData(stim, spikeTimes, groupIndex);
    
    modelParamsReal.regularize = 0;
    modelParamsReal.regularizeWeight = 1e6;
    
    modelParamsReal.debugL1 = 0;
    modelParamsReal.debugL2 = 0;
 
    tic;
    [modelParamsReal, gfd] = lnpsr2GradFD(modelParamsReal, datIdx);
    %gfd = zeros(1, nStimDelays*numStimChannels + nSpikeDelays);
    fdTime = toc;
    
    tic;
    [modelParamsReal, gFast] = lnpsr2Grad(modelParamsReal, datIdx);
    fastTime = toc;
    
    tic;
    [modelParamsReal, gLowMem] = lnpsr2GradLowMem(modelParamsReal, datIdx);
    lowMemTime = toc;
    
    
    gfd
    gFast
    gLowMem
    
    
    fastDiff = norm(gfd - gFast)
    lowMemDiff = norm(gfd - gLowMem)
    
    fdTime
    fastTime
    lowMemTime
    
    %{
    figure; hold on;
    plot(gfd, 'k-');
    plot(gFast, 'b-');
    plot(gLowMem, 'r-');
    legend('FD', 'Fast', 'LowMem');
    %}
    
    