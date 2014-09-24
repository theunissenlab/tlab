function test_glog()

    M = 0.25;
    B = 1.1;
    
    xvals = -5:0.1:5;
    yvals = glogistic(xvals, B, M);
    yvals = randn(size(yvals))*1e-1 + yvals;
    
    datIdx = 1:length(yvals);
    groupIndex = ones(size(yvals));
    
    strfData(xvals, yvals, groupIndex);
    
    gp = glogInit();
    xpnt = randn(1, 2);
    %gp.B = abs(xpnt(1)*2);
    %gp.M = abs(xpnt(2)*2);
    gp.B = 0.1;
    gp.M = 0;
    
    gp
    
    optOptions = trnSCG();
    optOptions.display = 1;
    
    [modelParamsTrained, optOptions] = strfOpt(gp, datIdx, optOptions);
    
    modelParamsTrained
    [modelParamsTrained, presp] = glogFwd(modelParamsTrained, datIdx);  
    
    
    figure; hold on;
    plot(xvals, yvals, 'k-');
    plot(xvals, presp, 'r-');
    legend('Actual', 'Predicted');
    
    