function test_gomp()

    a = 1;
    b = -5;
    c = -2;

    xvals = -5:0.1:5;
    yvals = gompertz(a, b, c, xvals);
    yvals = randn(size(yvals))*1e-2 + yvals;
    
    datIdx = 1:length(yvals);
    groupIndex = ones(size(yvals));
    
    strfData(xvals, yvals, groupIndex);
    
    gp = gompInit(a);
    xpnt = randn(1, 2);
    gp.b = -abs(xpnt(1)*2);
    gp.c = -abs(xpnt(2)*2);
    
    gp
    
    optOptions = trnSCG();
    optOptions.display = 1;
    
    [modelParamsTrained, optOptions] = strfOpt(gp, datIdx, optOptions);
    
    modelParamsTrained
    [modelParamsTrained, presp] = gompFwd(modelParamsTrained, datIdx);  
    
    
    figure; hold on;
    plot(xvals, yvals, 'k-');
    plot(xvals, presp, 'r-');
    legend('Actual', 'Predicted');
    
    