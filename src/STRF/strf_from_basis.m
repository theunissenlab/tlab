function strf = strf_from_basis(strfBasis)

    numChannels = strfBasis.numChannels;
    numDelays = strfBasis.numDelays;
    nChanFuncs = strfBasis.numChannelFunctions;
    nTempFuncs = strfBasis.numTemporalFunctions;
    channelSpacing = strfBasis.channelSpacing;
    temporalSpacing = strfBasis.temporalSpacing;
    
    weights = strfBasis.weights;

    strf = zeros(numChannels, numDelays);
    [xvals, yvals] = meshgrid(0:(numDelays-1), 0:(numChannels-1));
    
    for k = 1:nChanFuncs        
        for j = 1:nTempFuncs
            
            y0 = channelSpacing*(k-1) + 1;
            x0 = temporalSpacing*(j-1) + 1;
            w = weights(k, j);
            %fprintf('x0y0=(%d, %d), w=%f\n', x0, y0, w);
            params =  strfBasis.params(k, j, :);
            fval = strfBasis.basis.func(xvals, yvals, x0, y0, params);
            strf = strf + w*fval;
            
        end
    end
    
    