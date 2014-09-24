function test_strf_from_basis()

    numChannels = 60;
    numDelays = 75;
    
    channelSpacing = 3;
    temporalSpacing = 3;
    
    
    %% gaussian basis
    strfBasis = create_strf_basis(numChannels, numDelays, channelSpacing, temporalSpacing)    
    strfBasis.weights(5:10, [1:2:11]) = -1;
    strf = strf_from_basis(strfBasis);
    
    figure; hold on;
    imagesc(strf);
    axis tight;
    title('Gaussian2D');
    
    %% gabor basis
    gaborBasis = basis_gabor2d();
    strfBasis = create_strf_basis(numChannels, numDelays, channelSpacing, temporalSpacing, gaborBasis)
    strfBasis.weights(5:10, [1:2:11]) = -1;
    strf = strf_from_basis(strfBasis);
    
    figure; hold on;
    imagesc(strf);
    axis tight;
    title('Gabor2d');