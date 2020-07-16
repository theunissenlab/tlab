function test_glm()
    
    N = 20000;
    nChannels = 2;
    strfLen = 3;
    delays = 0:(strfLen-1);
    dispersion = 1;
    
    numTrials = 20;
    %family = init_glm_family_gaussian();
    %family = init_glm_family_binomial(numTrials);
    family = init_glm_family_poisson();
    
    %outputNL = @(x) identity_outputnl(x);
    %outputNL = @(x) sigmoid_outputnl(x);
    outputNL = @(x) exp_outputnl(x);
    
    datIdx = 1:N;
    
    
    %% initialize model
    modelParamsReal = glmInit(nChannels, delays, family, outputNL, dispersion);
    modelParamsReal.w1 = randn(nChannels, length(delays));
    modelParamsReal.b1 = randn();
    
    %% create stimulus
    stim = randn(N, nChannels);
    groupIndex = ones(1, N);
    resp = zeros(1, N);
    
    %% create response
    global globDat;
    strfData(stim, resp, groupIndex);
    
    [modelParamsReal, resp, linResp] = glmFwd(modelParamsReal, datIdx);
    strfData(stim, resp, groupIndex);    
    
    %% create model to fit
    modelParams = glmInit(nChannels, delays, family, outputNL, dispersion);
    %modelParams.w1 = randn(nChannels, length(delays));
        
    %% fit model
    optOptions = trnThreshGradDescLS();    
    optOptions.threshold = 0;
    optOptions.earlyStop = 0;
    optOptions.display = 1;
    optOptions.maxIter = 500;    
    optOptions.maxStep = 5;
    optOptions.minStep = 1e-7;
    optOptions.maxLSIter = 50;
    optOptions.stepSize = 1e-4;
    optOptions.stepConvWindow = 8;
    optOptions.stepConv = 5e-6;
    
    %% run direct fit
    [modelParamsTrained, optOptions] = strfOpt(modelParams, datIdx, optOptions);
    
    realWeights = modelParamsReal.w1;
    trainedWeights = modelParamsTrained.w1;
    wdiff = realWeights - trainedWeights;
    
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
    
    
    
end


function [y, dy] = exp_outputnl(x)
    y = exp(x);
    if nargout > 1
        dy = exp(x);
    end
end

function [y, dy] = sigmoid_outputnl(x)
    y = 1 ./ (1 + exp(-x));
    if nargout > 1
        dy = y .* (1 - y);
    end
end

function [y, dy] = identity_outputnl(x)
    y = x;
    if nargout > 1
        dy = 1;
    end
end