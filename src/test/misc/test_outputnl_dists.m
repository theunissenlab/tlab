function test_outputnl_dists()

    %% intialize a linear-nonlinear poisson model
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
    nGroups = 10;
    groupSize = nSamplePoints / nGroups;
    groupIndex = [];
    for k = 1:nGroups
        groupIndex = [groupIndex ones(1, groupSize)*k];
    end
    
    strfData(stim, resp, groupIndex);
    
    
    %% simulate spike trials to generate training data
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
    

    %% fit output nonlinearity
    nlinfoSpline = estimate_outputnl_from_spline(linResp, psth, numTrials, groupIndex);
    nlinfoDists = estimate_outputnl_from_dists(linResp, psth, numTrials, groupIndex);
    
    %% compute responses from p(spike|x) NL
    respFromP_err = fnval(nlinfoDists.err.outputNL, linResp);    
    respFromP_kl = fnval(nlinfoDists.kl.outputNL, linResp);

    
    %% compute responses from cubic spline approx
    respFromCs_err = fnval(nlinfoSpline.err.outputNL, linResp);
    respFromCs_kl = fnval(nlinfoSpline.kl.outputNL, linResp);
    
    
    %% preprocess model responses for info calculations
    infoFreqCutoff = 90;
    infoWindowSize = 0.250;
    respParams = struct;
    respParams.units = 's';
    respParams.split = 1;
    respParams.stimLengths = [nSamplePoints / sampleRate];
    [wholeResponse, groupIndex, respInfo, respParams] = preprocSpikeResponses(allSpikeTrials, respParams);
   
    
    %% compute performances
    [cBound, cP_err] = compute_coherence_full(respFromP_err, modelResponse, wholeResponse(1, :), wholeResponse(2, :), ...
                                              sampleRate, numTrials, infoFreqCutoff, infoWindowSize);
    pPerfRatio_err = cP_err.info / cBound.info
    
    [cBound, cP_kl] = compute_coherence_full(respFromP_kl, modelResponse, wholeResponse(1, :), wholeResponse(2, :), ...
                                              sampleRate, numTrials, infoFreqCutoff, infoWindowSize);
    pPerfRatio_kl = cP_kl.info / cBound.info    
    
    
    [cBound, cCs_err] = compute_coherence_full(respFromCs_err, modelResponse, wholeResponse(1, :), wholeResponse(2, :), ...
                                              sampleRate, numTrials, infoFreqCutoff, infoWindowSize);
    csPerfRatio_err = cCs_err.info / cBound.info
    
    
    [cBound, cCs_kl] = compute_coherence_full(respFromCs_err, modelResponse, wholeResponse(1, :), wholeResponse(2, :), ...
                                              sampleRate, numTrials, infoFreqCutoff, infoWindowSize);
    csPerfRatio_kl = cCs_kl.info / cBound.info
                                          
    
    %% plot simulation
    figure; hold on;
    subplot(3, 1, 1); hold on;
    plot(stim, 'k-');
    subplot(3, 1, 2); hold on;
    plot(fullResponse.stimCurrent, 'g-');
    subplot(3, 1, 3); hold on;
    plot(modelResponse, 'k-');
    
    %% plot output NL information for error objective function
    figure('Name', 'Lsq Error Obj Func');
    subplot(3, 1, 1); hold on;
    plot(nlinfoDists.err.x, nlinfoDists.err.px, 'k-'); axis tight;
    plot(nlinfoDists.err.x, nlinfoDists.err.px_smooth, 'k--'); axis tight;    
    plot(nlinfoDists.err.x, nlinfoDists.err.pxspike, 'r-'); axis tight;
    plot(nlinfoDists.err.x, nlinfoDists.err.pxspike_smooth, 'r--'); axis tight;    
    plot(nlinfoDists.err.x, nlinfoDists.err.pxnospike, 'b-'); axis tight;    
    plot(nlinfoDists.err.x, nlinfoDists.err.pxnospike_smooth, 'b--'); axis tight;    
    legend('p(x)', 'cs', 'p(x|spike)', 'cs', 'p(x|nospike)', 'cs');
    
    subplot(3, 1, 2); hold on;
    plot(nlinfoDists.err.x, nlinfoDists.err.pspikex, 'r-'); axis tight;
    plot(nlinfoDists.err.x, nlinfoDists.err.pspikex_smooth, 'r--'); axis tight;
    plot(nlinfoDists.err.x, nlinfoDists.err.pnospikex, 'b-'); axis tight;
    plot(nlinfoDists.err.x, nlinfoDists.err.pnospikex_smooth, 'b--'); axis tight;
    legend('p(spike|x)', 'cs', 'p(nospike|x)', 'cs');
    
    subplot(3, 1, 3); hold on;
    plot(nlinfoDists.err.x, fnval(nlinfoSpline.err.outputNL, nlinfoDists.err.x));
    axis tight;
    title('spline');
    
    
    %% plot output NL information for KL distance objective function
    figure('Name', 'KL Dist Obj Func');
    subplot(3, 1, 1); hold on;
    plot(nlinfoDists.kl.x, nlinfoDists.kl.px, 'k-'); axis tight;
    plot(nlinfoDists.kl.x, nlinfoDists.kl.px_smooth, 'k--'); axis tight;    
    plot(nlinfoDists.kl.x, nlinfoDists.kl.pxspike, 'r-'); axis tight;
    plot(nlinfoDists.kl.x, nlinfoDists.kl.pxspike_smooth, 'r--'); axis tight;    
    plot(nlinfoDists.kl.x, nlinfoDists.kl.pxnospike, 'b-'); axis tight;    
    plot(nlinfoDists.kl.x, nlinfoDists.kl.pxnospike_smooth, 'b--'); axis tight;    
    legend('p(x)', 'cs', 'p(x|spike)', 'cs', 'p(x|nospike)', 'cs');
    
    subplot(3, 1, 2); hold on;
    plot(nlinfoDists.kl.x, nlinfoDists.kl.pspikex, 'r-'); axis tight;
    plot(nlinfoDists.kl.x, nlinfoDists.kl.pspikex_smooth, 'r--'); axis tight;
    plot(nlinfoDists.kl.x, nlinfoDists.kl.pnospikex, 'b-'); axis tight;
    plot(nlinfoDists.kl.x, nlinfoDists.kl.pnospikex_smooth, 'b--'); axis tight;
    legend('p(spike|x)', 'cs', 'p(nospike|x)', 'cs');
    
    subplot(3, 1, 3); hold on;
    plot(nlinfoDists.kl.x, fnval(nlinfoSpline.kl.outputNL, nlinfoDists.kl.x));
    axis tight;
    title('spline');
    
    
    %% plot nonlinear responses
    figure; hold on;    
    subplot(4, 1, 1); hold on;
    plot(modelResponse, 'k-');
    plot(respFromP_err, 'b-');
    axis tight;
    title('NL from p(spike|x) (Lsq Err)');
    
    subplot(4, 1, 2); hold on;
    plot(modelResponse, 'k-');
    plot(respFromP_kl, 'b-');
    axis tight;
    title('NL from p(spike|x) (KL dist)');
    
    subplot(4, 1, 3); hold on;
    plot(modelResponse, 'k-');
    plot(respFromCs_err, 'b-');
    axis tight;
    title('NL from spline (Lsq err)');
    
    subplot(4, 1, 4); hold on;
    plot(modelResponse, 'k-');
    plot(respFromCs_kl, 'b-');
    axis tight;
    title('NL from spline (KL dist)');
    
    
    figure('Name', 'Talk Fig');
    subplot(2, 1, 1); hold on;
    plot(nlinfoDists.err.x, nlinfoDists.err.px, 'k-'); axis tight;    
    plot(nlinfoDists.err.x, nlinfoDists.err.pxspike, 'r-'); axis tight;    
    plot(nlinfoDists.err.x, nlinfoDists.err.pxnospike, 'b-'); axis tight;        
    legend('p(x)', 'p(x|spike)', 'p(x|nospike)');
    
    subplot(2, 1, 2); hold on;
    plot(nlinfoDists.err.x, nlinfoDists.err.pspikex, 'r-'); axis tight;    
    plot(nlinfoDists.err.x, nlinfoDists.err.pnospikex, 'b-'); axis tight;    
    legend('p(spike|x)', 'p(nospike|x)');
    
    
    
    
    