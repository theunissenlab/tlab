function test_estimate_spike_dists()

    %% intialize a linear-nonlinear poisson model to generate test data
    sampleRate = 1000;
    numStimChannels = 1;
    delays = 0:2;
    modelParams = glppInit(numStimChannels, delays, 'exponential', [], sampleRate, 0);
    modelParams.w1 = [1.8 1.5 -0.2];
    %modelParams.w1 = [2.8 1.5 0.2];
    modelParams.b1 = 0.001;
    
    nSamplePoints = 10000;
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
    
    
    %% estimate spike-triggered distributions
    distinfo = estimate_spike_dists(linResp, psth, numTrials);
    
    
    %% plot distributions and output nonlinearities
    figure; hold on;
    
    subplot(2, 1, 1); hold on;
    plot(distinfo.x, distinfo.px, 'k-'); axis tight;
    plot(distinfo.x, distinfo.px_smooth, 'k--'); axis tight;    
    plot(distinfo.x, distinfo.pxspike, 'r-'); axis tight;
    plot(distinfo.x, distinfo.pxspike_smooth, 'r--'); axis tight;    
    plot(distinfo.x, distinfo.pxnospike, 'b-'); axis tight;    
    plot(distinfo.x, distinfo.pxnospike_smooth, 'b--'); axis tight;    
    legend('p(x)', 'cs', 'p(x|spike)', 'cs', 'p(x|nospike)', 'cs');
    
    subplot(2, 1, 2); hold on;
    plot(distinfo.x, distinfo.pspikex, 'r-'); axis tight;
    plot(distinfo.x, distinfo.pspikex_smooth, 'r--'); axis tight;
    plot(distinfo.x, distinfo.pnospikex, 'b-'); axis tight;
    plot(distinfo.x, distinfo.pnospikex_smooth, 'b--'); axis tight;
    legend('p(spike|x)', 'cs', 'p(nospike|x)', 'cs');
    