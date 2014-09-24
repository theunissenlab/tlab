function strfBasis = fit_strf_basis(strf, strfBasis)

    objFunc = @(x) strfBasisErrorFunc(x, strfBasis, strf);
    
    opts = optimset('GradObj', 'on');
    
    w0 = randn(size(strfBasis.weights));
    
    wOpt = fminunc(objFunc, w0, opts);
    
    strfBasis.weights = reshape(wOpt, strfBasis.numChannelFunctions, strfBasis.numTemporalFunctions);
    