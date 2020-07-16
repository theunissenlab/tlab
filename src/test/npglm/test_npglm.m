function test_npglm()
    
    N = 20000;
    nChannels = 1;
    strfLen = 3;
    delays = 0:(strfLen-1);
    dispersion = 1;
    
    numTrials = 20;
    family = init_glm_family_gaussian();
    %family = init_glm_family_binomial(numTrials);
    %family = init_glm_family_poisson();
    
    datIdx = 1:N;
    
    
    %% initialize model
    ndim = 1;
    nb = 3;
    ctrs = linspace(-5, 5, nb)';
    
    modelParamsReal = npglmInit(nChannels, delays, family, dispersion, @exp_bfunc, ctrs);
    modelParamsReal.w1 = randn(nChannels, length(delays));
    modelParamsReal.w2 = randn(1, nb);
    modelParamsReal.b1 = randn();
    
    %% create stimulus
    stim = randn(N, nChannels);
    groupIndex = ones(1, N);
    resp = zeros(1, N);
    
    %% create response
    global globDat;
    strfData(stim, resp, groupIndex);
    
    [modelParamsReal, resp] = npglmFwd(modelParamsReal, datIdx);
    strfData(stim, resp, groupIndex);    
    
    %% create model to fit
    modelParams = npglmInit(nChannels, delays, family, dispersion, @exp_bfunc, ctrs);
    
    wgain = 1e-1;
    %{
    modelParams.w1 = modelParamsReal.w1 + randn(nChannels, length(delays))*wgain;
    modelParams.w2 = modelParamsReal.w2 + randn(1, nb)*wgain;
    modelParams.b1 = modelParamsReal.b1 + randn()*wgain;
    %}
    modelParams.w1 = randn(nChannels, length(delays))*wgain;
    modelParams.w2 = randn(1, nb)*wgain;
    modelParams.b1 = randn()*wgain;
    
    
    %% fit model    
    
    optOptions = trnThreshGradDescLS();    
    optOptions.threshold = 0;
    optOptions.earlyStop = 0;
    optOptions.display = 1;
    optOptions.maxIter = 400;
    optOptions.maxStep = 5;
    optOptions.minStep = 1e-7;
    optOptions.maxLSIter = 50;
    optOptions.stepSize = 1e-4;
    optOptions.stepConvWindow = 8;
    optOptions.stepConv = 1e-8;
    
    %{
    optOptions = trnSCG();    
    optOptions.maxIter = 500;
    optOptions.display = 1;
    %}
    
    %% train
    [modelParamsTrained, optOptions] = strfOpt(modelParams, datIdx, optOptions);
    
    realWeights = modelParamsReal.w1
    trainedWeights = modelParamsTrained.w1
    wdiff = realWeights - trainedWeights
    
    realBasisWeights = modelParamsReal.w2
    trainedBasisWeights = modelParamsTrained.w2
    wdiff2 = realBasisWeights - trainedBasisWeights
    
    %{
    figure; hold on;
    subplot(3, 1, 1); hold on;
    imagesc(realWeights); axis tight; colorbar;
    title(sprintf('bias=%f', modelParamsReal.b1));
    subplot(3, 1, 2); hold on;
    imagesc(trainedWeights); axis tight; colorbar;
    title(sprintf('bias=%f', modelParamsTrained.b1));
    subplot(3, 1, 3); hold on;
    imagesc(wdiff); axis tight; colorbar;
    title(sprintf('norm diff=%f', norm(wdiff)));
    %}
end


function [y, dy] = exp_bfunc(r)

    y = exp(-(r.^2));    
    dy = -2 .* r .* exp(-(r.^2));
end
