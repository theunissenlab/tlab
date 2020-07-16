function test_glm_grad()
    
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
    modelParams = glmInit(nChannels, delays, family, outputNL, dispersion);
    modelParams.w1 = randn(nChannels, length(delays));
    modelParams.b1 = 0;
    
    %% create stimulus
    stim = randn(N, nChannels);
    groupIndex = ones(1, N);
    resp = zeros(1, N);
    
    %% create response
    global globDat;
    strfData(stim, resp, groupIndex);
    
    [modelParams, resp, linResp] = glmFwd(modelParams, datIdx);
    strfData(stim, resp, groupIndex);    
    
    %% perturb weights of model
    modelParams.w1 = modelParams.w1 + randn(nChannels, length(delays));
    modelParams.b1 = modelParams.b1 + randn();
    
    %% compute gradient
    [modelParams, g] = glmGrad(modelParams, datIdx);

    %% compute FD approximation to gradient
    [modelParams, gfd] = glmGradFD(modelParams, datIdx);
    
    g
    gfd
    
    
    gdiff = g - gfd;
    
    figure; hold on;
    subplot(3, 1, 1); hold on;
    imagesc(reshape(g(1:end-1), nChannels, strfLen)); axis tight; colorbar;    
    subplot(3, 1, 2); hold on;
    imagesc(reshape(gfd(1:end-1), nChannels, strfLen)); axis tight; colorbar;
    subplot(3, 1, 3); hold on;
    imagesc(reshape(gdiff(1:end-1), nChannels, strfLen)); axis tight; colorbar;
    title(sprintf('norm diff=%f', norm(gdiff)));
    
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
