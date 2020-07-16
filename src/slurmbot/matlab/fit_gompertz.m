function gompParams = fit_gompertz(xvals, yvals, a, b0, c0, trainingIndex)

    groupIndex = ones(size(yvals));

    strfData(rv(xvals), rv(yvals), groupIndex);
        
    gp = gompInit(a);
    gp.b = b0;
    gp.c = c0;
    
    optOptions = trnSCG();
    optOptions.display = 1;
    
    [modelParamsTrained, optOptions] = strfOpt(gp, trainingIndex, optOptions);
    
    gompParams = modelParamsTrained;