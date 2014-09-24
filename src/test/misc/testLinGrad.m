function testLinGrad()

    duration = 1;
    sampleRate = 1000;
    nDataPoints = round(duration*sampleRate);

    nStimDelays = 40;
    numStimChannels = 2;
    stimDelays = 0:(nStimDelays-1);
    
    nlType = 'linear';
    %nlType = 'exponential';
    %nlType = 'softmax';
    %nlType = 'logistic';

    %% initialize model
    modelParamsReal = linInit(numStimChannels, stimDelays, nlType);
    %modelParamsReal.freqDomain = 1;
    
    modelParamsReal.w1 = randn(numStimChannels, nStimDelays);
    modelParamsReal.b1 = 0.5;
    
    %% create test data
    stim = randn(nDataPoints, numStimChannels);

    % set up global structure
    groupIndex = ones(nDataPoints, 1);
    global globDat;
    strfData(stim, zeros(1, nDataPoints), groupIndex);

    datIdx = 1:nDataPoints;
    [modelParamsReal, modelResponseReal] = linFwd(modelParamsReal, datIdx);
    
    strfData(stim, modelResponseReal, groupIndex);
    
    %randomize strf to get a nonzero error
    modelParamsReal.w1 = randn(numStimChannels, nStimDelays);
    
    tic
    gfd = linGradFD(modelParamsReal, datIdx);
    %gfd = zeros(1, nStimDelays*numStimChannels + nSpikeDelays);
    fdTime = toc;
    
    tic;
    [modelParamsReal, gReal] = linGrad(modelParamsReal, datIdx);
    realTime = toc;
    
    gfd
    gReal
    
    gDiff = norm(gfd - gReal)
    
    fdTime
    realTime
    
    figure; hold on;
    plot(gfd-gReal, 'r-');
    title('Difference in Gradients');
    