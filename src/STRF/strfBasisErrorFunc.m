function [err, grad] = strfBasisErrorFunc(weights, strfBasis, strf)
    
    strfBasis.weights = reshape(weights, strfBasis.numChannelFunctions, strfBasis.numTemporalFunctions);
        
    strfFromBasis = strf_from_basis(strfBasis);
    
    err = sum(sum((strf - strfFromBasis).^2));

    if nargout > 1

        numChannels = strfBasis.numChannels;
        numDelays = strfBasis.numDelays;
        numFuncs = strfBasis.numFunctions;
        nChanFuncs = strfBasis.numChannelFunctions;
        nTempFuncs = strfBasis.numTemporalFunctions;
        channelSpacing = strfBasis.channelSpacing;
        temporalSpacing = strfBasis.temporalSpacing;
        
        linStrf = strf(:);
        linStrfFromBasis = strfFromBasis(:);
        
        numTimePoints = length(linStrf);
        
        e = linStrf - linStrfFromBasis;        
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
        
        grad = -2*(phi' * e);
        
    end
    