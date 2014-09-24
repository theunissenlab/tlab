function test_lnpsr

    nDataPoints = 1000;

    nStimDelays = 3;
    numStimChannels = 2;
    stimDelays = 0:(nStimDelays-1);
        
    %nSpikeDelays = 9;
    %nSpikeDelays = 4;
    %spikeDelays = 0:(nSpikeDelays-1);
    nSpikeDelays = 0;
    spikeDelays = [];
    
    
    nlType = 'exponential';
    
    sampleRate = 1000;
    
    numTrials = 20;
    %% initialize model
    
    modelParamsReal = lnpsrInit(numStimChannels, stimDelays, nlType, spikeDelays, sampleRate, numTrials);
    
    %w = [0.25 0.4 0.3; 0.1 0.2 0.5];
    w = [0.25 0.10 0.05; 0.15 0.05 0.01];
    %w = [2];
    b = 0.1;    
    
    modelParamsReal.w1 = w;
    modelParamsReal.b1 = b;
    
    %modelParamsReal.spikeResponseWeights = [-1.5 -1.5 -1 -0.5];    
    %modelParamsReal.spikeResponseWeights = [1.5 1 0.5 0.25 0 -0.5 -0.5 -0.25 -0.1];    
    modelParamsReal.spikeResponseWeights = [];
    
    %% create test data
    gain = 5;
    stim = gen_test_stim(nDataPoints, numStimChannels, sampleRate)*gain;
    
    %subtract off mean stimulus
    %stim = stim - mean(stim);
    
    % set up global structure
    groupIndex = ones(nDataPoints, 1);
    global globDat;
    strfData(stim, zeros(1, nDataPoints), groupIndex);
    
    datIdx = 1:nDataPoints;
    [modelParamsReal, modelResponseReal, fullResponse] = lnpsrFwd(modelParamsReal, datIdx, numTrials);
    modelResponseReal = modelResponseReal*numTrials;
    
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
            plot2svg_2d('~/temp/lnpsr-linresp.svg', h1);

            h2 = figure; hold on;
            spIndx = find(fullResponse.spikeTrials(k, :) > 0);
            %subplot(2, 1, 2); hold on;
            plot(t, fullResponse.nonlinearResponse(k, :), 'b-');
            plot(t(spIndx), ones(1, length(spIndx)), 'ro');
            title(sprintf('%d: Nonlinear Response', k));
            %axis([min(t) max(t) 0 5]);
            axis tight;
            plot2svg_2d('~/temp/lnpsr-nonlinresp.svg', h2);
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
    
    optOptions.display = 1;
        
    %% set up global structure with fake data
    groupIndex = ones(nDataPoints, 1);
    strfData(stim, modelResponseReal, groupIndex);
    
    %% Initialize LNP spike-response model
    nlType = 'exponential';
    modelParams = lnpsrInit(numStimChannels, stimDelays, nlType, spikeDelays, sampleRate, numTrials);
    modelParams.l2weight = 1;
        
    %% run strflab optimization
    [modelParamsTrained, optOptions] = strfOpt(modelParams, datIdx, optOptions);
       
    %% look at response
    [modelParamsTrained, predResp] = lnpsrFwd(modelParamsTrained, datIdx, numTrials);
    predResp = (predResp/max(predResp)) * max(modelResponseReal);
    
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
