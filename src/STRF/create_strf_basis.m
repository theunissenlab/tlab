function strfBasis = create_strf_basis(numChannels, numDelays, channelSpacing, temporalSpacing, basis)

    if nargin < 5
        basis = basis_gaussian2d();
    end
    
    nCWts = floor(numChannels / channelSpacing) + 1;
    nTWts = floor(numDelays / temporalSpacing) + 1;
    
    weights = ones(nCWts, nTWts);
    
    nChanFuncs = size(weights, 1);
    nTempFuncs = size(weights, 2);
    nFuncs = nChanFuncs*nTempFuncs;
    
    params = zeros(nChanFuncs, nTempFuncs, basis.numParams);
    
    strfBasis = struct;
    strfBasis.numChannels = numChannels;
    strfBasis.numDelays = numDelays;
    strfBasis.channelSpacing = channelSpacing;
    strfBasis.temporalSpacing = temporalSpacing;
    strfBasis.basis = basis;
    strfBasis.numChannelFunctions = nChanFuncs;
    strfBasis.numTemporalFunctions = nTempFuncs;
    strfBasis.numFunctions = nFuncs;
    strfBasis.weights = weights;    
    strfBasis.params = params;
    strfBasis.numTotalParams = nFuncs + nFuncs*basis.numParams;
    