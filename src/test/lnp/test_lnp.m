function test_lnp

    nDataPoints = 1000;

    nDelays = 1;
    nChannels = 1;
    delays = 0:(nDelays-1);
    modelParamsReal = lnpInit(nChannels, delays);
        
    %w = [0.25 0.4 0.3; 0.1 0.2 0.5];
    w = [0.25];
    
    modelParamsReal.w1 = w;
    modelParamsReal.b1 = 0.0;
    
    spikeTimes = zeros(nDataPoints, 1);
    spikeIndicies = [5 10 15 24 36 41 43 55 60 80 90];
    spikeIndicies = round((spikeIndicies / 100)*nDataPoints);
    spikeTimes(spikeIndicies) = 1;
    
    stim = randn(nDataPoints, nChannels);
    stim(spikeIndicies) = 15;
    
    groupIndex = ones(nDataPoints, 1);
    
    global globDat;
    strfData(stim, spikeTimes, groupIndex);

    datIdx = [];
    [modelParamsReal, modelResponseReal, modelResponseLinear] = lnpFwd(modelParamsReal, datIdx);
    
    %% Set options for optimization method
    %{
    optOptions = trnGradDescLineSearch;
    optOptions.display = 1;
    optOptions.coorDesc = 0;
    optOptions.earlyStop = 0;
    optOptions.stepSize = 1e-5;
    optOptions.nDispTruncate = 0;
    optOptions.maxIter = 10000;
    optOptions.lineSearch = 1;
    %}
    %optOptions = trnSCGNetlab;
    optOptions = trnSCG;
    
    optOptions.display = 1;
    
    
    %% Initialize linear model
    modelParams = lnpInit(size(stim,2), 0:(nDelays-1));
    modelParams.w1 = randn(nChannels, length(delays));
    modelParams.e2Weight = 1;

    %% set the output nonlinearity
    modelParams.outputNL = 'exponential';
        
    %% run strflab optimization
    [modelParamsTrained, optOptions] = strfOpt(modelParams, datIdx, optOptions);
       
    %% look at response
    [modelParamsTrained, predResp] = strfFwd(modelParamsTrained, datIdx);
    predResp = (predResp / max(predResp))*max(modelResponseReal);
    
    realStrf = modelParamsReal.w1
    predStrf = modelParamsTrained.w1
    strfDiff = norm(realStrf - predStrf)
  
    realw0 = modelParamsReal.b1
    predw0 = modelParamsTrained.b1
    
    figure; hold on;
    plot(modelResponseReal, 'k-', 'LineWidth', 3);
    plot(predResp, 'r-');
    legend('Real', 'Model');
    