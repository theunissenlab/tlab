function test_outputnl_joint()

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
    

    %% fit output nonlinearity using 3 different methods    
    [nlinfoJoint, nlparamsJoint] = estimate_outputnl_from_joint(linResp, psth, numTrials);    
    
    
    %% compute responses for joint NL
    
    jointModelParams = glppInit(numStimChannels, delays, 'custom', [], sampleRate, 0);
    jointModelParams.w1 = modelParams.w1;
    jointModelParams.b1 = modelParams.b1;
    jointModelParams.nlFunc = @(x) fnval(nlinfoJoint.outputNL, x)*sampleRate;    
    [jointModelParams, respFromJoint] = glppFwd(jointModelParams, datIdx, numTrials);
    %respFromJoint = zeros(1, length(datIdx));
    
    
    %% preprocess model responses for info calculations
    infoFreqCutoff = 90;
    infoWindowSize = 0.250;
    respParams = struct;
    respParams.units = 's';
    respParams.split = 1;
    respParams.stimLengths = [nSamplePoints / sampleRate];
    [wholeResponse, groupIndex, respInfo, respParams] = preprocSpikeResponses(allSpikeTrials, respParams);
   
    
    %% compute performances
    [cBound, cJoint] = compute_coherence_full(respFromJoint, modelResponse, wholeResponse(1, :), wholeResponse(2, :), ...
                                              sampleRate, numTrials, infoFreqCutoff, infoWindowSize);
    jointPerfRatio = cJoint.info / cBound.info
    
                                          
    
    %% plot simulation
    figure; hold on;
    subplot(3, 1, 1); hold on;
    plot(stim, 'k-');
    subplot(3, 1, 2); hold on;
    plot(fullResponse.stimCurrent, 'g-');
    subplot(3, 1, 3); hold on;
    plot(modelResponse, 'k-');
    
    %% plot distributions and output nonlinearities
    figure; hold on;
    
    subplot(4, 1, 1); hold on;
    plot(nlinfoDists.x, nlinfoJoint.px, 'k-'); axis tight;
    plot(nlinfoDists.x, nlinfoDists.pxspike, 'r-'); axis tight;
    plot(nlinfoDists.x, nlinfoDists.pxnospike, 'b-'); axis tight;    
    legend('p(x)', 'p(x|spike)', 'p(x|nospike)');
    
    subplot(4, 1, 3); hold on;
    plot(nlinfoJoint.x, exp(nlinfoJoint.x)/1000, 'k-');
    plot(nlinfoJoint.x, fnval(nlinfoJoint.outputNL, nlinfoJoint.x), 'r-');
    plot(nlinfoJoint.x, nlinfoJoint.rawNL, 'b-');
    axis tight;
    legend('Real', 'Pred', 'Raw');
    title('Output NL');
    
    subplot(4, 1, 4); hold on;
    imagesc(nlinfoJoint.x, nlinfoJoint.y, nlinfoJoint.pxy);
    %plot(nlinfoJoint.pxy_data(:, 1), nlinfoJoint.pxy_data(:, 2), 'w.');
    axis([min(nlinfoJoint.x) max(nlinfoJoint.x) min(nlinfoJoint.y) max(nlinfoJoint.y)]);
    
    
    %% plot nonlinear responses
    figure; hold on;    
    subplot(3, 1, 1); hold on;
    plot(modelResponse, 'k-');
    plot(respFromJoint, 'b-');
    axis tight;
    title('NL from Joint');
    
    subplot(3, 1, 2); hold on;
    plot(modelResponse, 'k-');
    plot(respFromP, 'b-');
    axis tight;
    title('NL from p(spike|x)');
    
    subplot(3, 1, 3); hold on;
    plot(modelResponse, 'k-');
    plot(respFromCs, 'b-');
    axis tight;
    title('NL from spline');
    %}