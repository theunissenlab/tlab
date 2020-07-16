function test_outputnl()

    %% intialize a linear-nonlinear poisson model to generate test data
    numStimChannels = 1;
    delays = 0:2;
    modelParams = glppInit(numStimChannels, delays, 'exponential', [], 1000, 0);
    modelParams.w1 = [1.8 1.5 -0.2];
    %modelParams.w1 = [2.8 1.5 0.2];
    modelParams.b1 = 0.001;
    
    stim = randn(30000, 1);
    resp = stim;
    groupIndex = ones(size(stim));
    strfData(stim, resp, groupIndex);
    
    numTrials = 50;
    datIdx = 1:length(groupIndex);
    [modelParams, modelResponse, fullResponse] = glppFwd(modelParams, datIdx, numTrials);
    
    figure; hold on;
    subplot(3, 1, 1); hold on;
    plot(stim, 'k-');
    subplot(3, 1, 2); hold on;
    plot(fullResponse.stimCurrent, 'g-');
    subplot(3, 1, 3); hold on;
    plot(modelResponse, 'k-');
    
    linResp = fullResponse.stimCurrent;
    psth = modelResponse;
    
    nlinfo = get_outputnl_info(linResp, psth, numTrials);
    
    [xvals, yvals, stds] = get_binned_outputnl(linResp, psth, 1);
    
    linStart = max(find(nlinfo.linRng <= min(xvals)));
    linEnd = min(find(nlinfo.linRng >= max(xvals)));
    linIndx = linStart:linEnd;
    
    figure; hold on;
    subplot(5, 1, 1); hold on;
    plot(nlinfo.linRng(linIndx), nlinfo.px(linIndx), 'k-'); axis tight;
    title(sprintf('p(x), bw=%0.4f', nlinfo.bwpx));    
    
    subplot(5, 1, 2); hold on;
    plot(nlinfo.linRng(linIndx), nlinfo.pxspike(linIndx), 'k-'); axis tight;
    title(sprintf('p(x|spike), bw=%0.4f', nlinfo.bwpxspike));
    
    subplot(5, 1, 3); hold on;
    plot(nlinfo.linRng(linIndx), nlinfo.pspikex(linIndx), 'k-'); axis tight;
    evals = exp(nlinfo.linRng(linIndx)) / 1000;
    plot(nlinfo.linRng(linIndx), evals, 'r-');
    %plot(nlinfo.linRng(linIndx), 1-exp(-evals), 'g-');
    plot(nlinfo.linRng(linIndx), nlinfo.outputNL(linIndx), 'b-');    
    plot(nlinfo.linpsth.linVals, nlinfo.pspikex_ds, 'g-');
    axis([min(nlinfo.linRng(linIndx)) max(nlinfo.linRng(linIndx)) 0 1]);
    title('p(spike|x)');
        
    subplot(5, 1, 4); hold on;
    plot(nlinfo.linpsth.linVals, nlinfo.linpsth.psthPreds, 'k-'); axis tight;
    plot(nlinfo.linpsth.linVals, nlinfo.linpsth.psthPredsCorr, 'g-'); axis tight;
    evals = exp(nlinfo.linpsth.linVals) / 1000;
    plot(nlinfo.linpsth.linVals, evals, 'r-'); axis tight;        
    axis([min(nlinfo.linpsth.linVals) max(nlinfo.linpsth.linVals) 0 max(nlinfo.linpsth.psthPreds)]);
    title('E[psth|x]');
    
    subplot(5, 1, 5); hold on;        
    imagesc(nlinfo.linpsth.linVals, nlinfo.linpsth.psthVals, nlinfo.linpsth.joint); axis tight;    
    