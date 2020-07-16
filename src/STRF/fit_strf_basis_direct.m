function strfBasis = fit_strf_basis_direct(strf, strfBasis)

    numChannels = strfBasis.numChannels;
    numDelays = strfBasis.numDelays;
    numFuncs = strfBasis.numFunctions;
    nChanFuncs = strfBasis.numChannelFunctions;
    nTempFuncs = strfBasis.numTemporalFunctions;
    channelSpacing = strfBasis.channelSpacing;
    temporalSpacing = strfBasis.temporalSpacing;

    linStrf = strf(:);

    numTimePoints = length(linStrf);

    phi = zeros(numTimePoints, numFuncs);

    [xvals, yvals] = meshgrid(0:(numDelays-1), 0:(numChannels-1));

    colIndx = 1;        
    for j = 1:nTempFuncs
        for k = 1:nChanFuncs        
            y0 = channelSpacing*(k-1) + 1;
            x0 = temporalSpacing*(j-1) + 1;
            params =  strfBasis.params(k, j, :);
            fval = strfBasis.basis.func(xvals, yvals, x0, y0, params);

            phi(:, colIndx) = fval(:);
            colIndx = colIndx + 1;
        end
    end

    strfLin = strf(:);
    
    phiSquared = phi' * phi;
    
    wOpt = phiSquared \ (phi' * strfLin);
    
    strfBasis.weights = reshape(wOpt, strfBasis.numChannelFunctions, strfBasis.numTemporalFunctions);
    