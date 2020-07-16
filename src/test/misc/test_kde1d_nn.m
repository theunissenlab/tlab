function test_kde1d_nn()

    %% intialize a linear-nonlinear poisson model to generate test data
    sampleRate = 1000;
    numStimChannels = 1;
    delays = 0:2;
    modelParams = glppInit(numStimChannels, delays, 'exponential', [], sampleRate, 0);
    modelParams.w1 = [1.8 1.5 -0.2];
    %modelParams.w1 = [2.8 1.5 0.2];
    modelParams.b1 = 0.001;
    
    nSamplePoints = 1000;
    stim = randn(nSamplePoints, 1);
    resp = stim;
    groupIndex = ones(size(stim));
    strfData(stim, resp, groupIndex);
    
    
    %% simulate spike trials
    numTrials = 20;
    datIdx = 1:length(groupIndex);
    [modelParams, modelResponse, fullResponse] = glppFwd(modelParams, datIdx, numTrials);
    
    linResp = fullResponse.stimCurrent;
    psth = modelResponse;
    
    spikeTrials = cell(numTrials, 1);
    for k = 1:numTrials       
        stimes = (find(fullResponse.spikeTrials(k, :))-1) / sampleRate;
        spikeTrials{k} = stimes;        
    end    
    allSpikeTrials = {spikeTrials};
    
    
    
    %% fit p(x)
    x = linResp;
    xmin = min(linResp);
    xmax = max(linResp);
    xinc = (xmax - xmin) / 250;    
    xgrid = xmin:xinc:xmax;    
    px = kde1d_nn(x, xgrid)
    
    figure;
    plot(xgrid, px, 'k-'); axis tight;
    