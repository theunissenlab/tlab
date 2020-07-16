function test_leglm_grad()
    
    N = 20000;
    nChannels = 2;
    strfLen = 3;
    delays = 0:(strfLen-1);
    dispersion = 1;
    
    numTrials = 20;
    %family = init_glm_family_gaussian();
    family = init_glm_family_binomial(numTrials);
    %family = init_glm_family_poisson();
        
    datIdx = 1:N;
    
    %% intialize rbf model
    
    %% initialize leglm model
    modelParams = leglmInit(nChannels, delays, family, dispersion);
    modelParams.w1 = randn(nChannels, length(delays));
    modelParams.b1 = randn();
    modelParams.m = rand()*4;
        
    %% create stimulus
    stim = randn(N, nChannels);
    groupIndex = ones(1, N);
    resp = zeros(1, N);
    
    %% create response
    global globDat;
    strfData(stim, resp, groupIndex);
    
    [modelParams, resp] = leglmFwd(modelParams, datIdx);    
    strfData(stim, resp, groupIndex);    
    
    %% perturb weights of model    
    modelParams.w1 = modelParams.w1 + randn(nChannels, length(delays));
    modelParams.b1 = modelParams.b1 + randn();
    modelParams.m = modelParams.m + rand();
        
    %% compute gradient
    [modelParams, g] = leglmGrad(modelParams, datIdx);

    %% compute FD approximation to gradient
    [modelParams, gfd] = leglmGradFD(modelParams, datIdx);
    
    g
    gfd
    
    gdiff = g - gfd
    
    ngdiff_lin = norm(g(end-1) - gfd(end-1))
    
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
