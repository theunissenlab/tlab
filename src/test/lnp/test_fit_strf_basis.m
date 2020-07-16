function test_fit_strf_basis()

    %% generate artificial STRF basis
    numChannels = 60;
    numDelays = 75;
    
    channelSpacing = 3;
    temporalSpacing = 3;
    
    % gabor basis
    gaborBasis = basis_gabor2d();
    strfBasis = create_strf_basis(numChannels, numDelays, channelSpacing, temporalSpacing, gaborBasis)
    strfBasis.weights = randn(size(strfBasis.weights));
    strf = strf_from_basis(strfBasis);
    
    %% test gradient
    %{
    w0 = randn(size(strfBasis.weights));
    w0Lin = w0(:);
    [err, grad] = strfBasisErrorFunc(w0Lin, strfBasis, strf);
    %finite diff approx
    deps = 1e-8;
    gradFd = zeros(size(grad));
    
    for k = 1:length(grad)
     
        w0Back = w0Lin;
        w0Fwd = w0Lin;
        w0Back(k) = w0Back(k) - deps;
        w0Fwd(k) = w0Fwd(k) + deps;
        
        errFwd = strfBasisErrorFunc(w0Fwd, strfBasis, strf);
        errBack = strfBasisErrorFunc(w0Back, strfBasis, strf);
        
        gradFd(k) = (errFwd - errBack) / (2*deps);
    end
    
    figure; hold on;
    plot(grad, 'k-');
    plot(gradFd, 'r-');
    axis tight;
    legend('Analytical', 'FD');
    ndiff = norm(grad - gradFd);
    title(sprintf('ndiff=%f', ndiff));
    %}
    
    %% test fitting
    %{
    strfBasisFit = create_strf_basis(numChannels, numDelays, channelSpacing, temporalSpacing, gaborBasis);
    %strfBasisFit = fit_strf_basis(strf, strfBasisFit);
    strfBasisFit = fit_strf_basis_direct(strf, strfBasisFit);
        
    strfFit = strf_from_basis(strfBasisFit);
    sdiff = norm(strfFit(:) - strf(:))
    
    figure; hold on;
    subplot(1, 2, 1);
    imagesc(strf);
    axis tight;
    
    subplot(1, 2, 2);
    imagesc(strfFit);
    axis tight;
    %}
    
    %% test fitting 2
    svars = load('~/berkeley/tlab/trunk/src/test/lnp/sample_strf.mat');
    strf = svars.strf;
    clear svars;
    
    sz = size(strf);
    
    strfBasisFit = create_strf_basis(sz(1), sz(2), channelSpacing, temporalSpacing, gaborBasis);
    strfBasisFit = fit_strf_basis_direct(strf, strfBasisFit);
        
    strfFit = strf_from_basis(strfBasisFit);
    sdiff = norm(strfFit(:) - strf(:))
    
    wLin = strfBasisFit.weights(:);
    wMean = mean(wLin);
    wStd = std(wLin);
    
    wLinFilt = wLin;
    wLinFilt(abs(wLin) < wStd) = 0;
    wFilt = reshape(wLinFilt, strfBasisFit.numChannelFunctions, strfBasisFit.numTemporalFunctions);
    
    strfBasisFilt = strfBasisFit;
    strfBasisFilt.weights = wFilt;
    strfFilt = strf_from_basis(strfBasisFilt);
    
    numFuncs = sum(wLinFilt > 0);
    
    %% plot
    figure; hold on;
    subplot(2, 2, 1);
    imagesc(strf);
    axis tight;
    title('Real STRF');
    
    subplot(2, 2, 2);
    imagesc(strfFit);
    axis tight;
    title('Basis STRF');
    
    subplot(2, 2, 3);
    hist(wLin, 35);
    title(sprintf('Weights: mean=%f, std=%f', wMean, wStd));
        
    subplot(2, 2, 4);
    imagesc(strfFilt);
    axis tight;
    title(sprintf('Filtered STRF: # funcs=%d', numFuncs));
       