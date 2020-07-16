function test_npglm_grad()
    
    N = 20000;
    nChannels = 2;
    strfLen = 3;
    delays = 0:(strfLen-1);
    dispersion = 1;
    
    numTrials = 20;
    family = init_glm_family_gaussian();
    %family = init_glm_family_binomial(numTrials);
    %family = init_glm_family_poisson();
        
    datIdx = 1:N;
    
    %% intialize rbf model
    ndim = 1;
    nb = 5;
    ctrs = linspace(-5, 5, nb)';
    
    %% initialize npglm model
    modelParams = npglmInit(nChannels, delays, family, dispersion, @exp_bfunc, ctrs);
    %modelParams = npglmInit(nChannels, delays, family, dispersion, @identity_bfunc, ctrs);
    modelParams.w1 = randn(nChannels, length(delays));
    modelParams.w2 = randn(1, nb);
    modelParams.b1 = randn();
    
    %% create stimulus
    stim = randn(N, nChannels);
    groupIndex = ones(1, N);
    resp = zeros(1, N);
    
    %% create response
    global globDat;
    strfData(stim, resp, groupIndex);
    
    [modelParams, resp] = npglmFwd(modelParams, datIdx);
    strfData(stim, resp, groupIndex);    
    
    %% perturb weights of model
    modelParams.w1 = modelParams.w1 + randn(nChannels, length(delays));
    modelParams.b1 = modelParams.b1 + randn();
    modelParams.w2 = modelParams.w2 + randn(1, nb);
    
    %% compute gradient
    [modelParams, g] = npglmGrad(modelParams, datIdx);

    %% compute FD approximation to gradient
    [modelParams, gfd] = npglmGradFD(modelParams, datIdx);
    
    g
    gfd
    
    gdiff = g - gfd
    
    ngdiff = norm(gdiff)
    
    %{
    figure; hold on;
    subplot(3, 1, 1); hold on;
    imagesc(reshape(g(1:end-1), nChannels, strfLen)); axis tight; colorbar;    
    subplot(3, 1, 2); hold on;
    imagesc(reshape(gfd(1:end-1), nChannels, strfLen)); axis tight; colorbar;
    subplot(3, 1, 3); hold on;
    imagesc(reshape(gdiff(1:end-1), nChannels, strfLen)); axis tight; colorbar;
    title(sprintf('norm diff=%f', norm(gdiff)));
    %}
    
end


function [y, dy] = exp_bfunc(r)

    y = exp(-(r.^2));    
    dy = -2 .* r .* exp(-(r.^2));
end

function [y, dy] = identity_bfunc(r)
    y = r;
    dy = ones(size(r));
end
