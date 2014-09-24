function fit_outputnl_once(outputFileName)

    vars = load(outputFileName);
    
    modelParamsTrained = vars.modelParamsTrained;
    preprocData = vars.preprocData;
    perfData = vars.perfData;
    trainingIndex = vars.trainingIndex;
    validationIndex = vars.validationIndex;
    
    clear vars;
    
    groupIndex = preprocData.groupIndex;
    wholeStim = preprocData.stim;
    wholeResponse = preprocData.resp;
    
    strfData(wholeStim, wholeResponse, groupIndex);
        
    realResp = perfData.wholeResp(trainingIndex);
    
    [modelParamsTrained, respTrain] = strfFwd(modelParamsTrained, trainingIndex);
    respTrain(isnan(respTrain)) = 0;
    
    [xvals, yvals, stds] = get_binned_outputnl(respTrain, realResp);
    
    size(yvals)
    size(xvals)
    
    slopes = diff(yvals) ./ diff(xvals)
    
    %{
    B = 1.75;
    M = 1.67;
    glogy = glogistic(xvals, B, M);
        
    gb = -7;
    gc = -1.5;
    gompy = gompertz(1, gb, gc, xvals);    
    
    figure; hold on;
    plot(xvals, yvals, 'ko-');
    plot(xvals, glogy, 'r-');
    plot(xvals, gompy, 'b-');
    title(sprintf('B=%f | M=%f || gb=%f | gc=%f', B, M, gb, gc));
    %}
    
    [B, M] = fit_glog_lsq(xvals, yvals, abs(randn()*2), 0);
    glogy = glogistic(xvals, B, M);
    
    agomp = 1;
    [gb, gc] = fit_gompertz_fmin(xvals, yvals, agomp, -abs(randn()), -abs(randn()));
    gompy = gompertz(agomp, gb, gc, xvals);
    
    
    figure; hold on;
    plot(xvals, yvals, 'ko-');
    plot(xvals, glogy, 'r-');
    plot(xvals, gompy, 'b-');
    legend('Orig', 'glog', 'gomp');
    title(sprintf('B=%f | M=%f || gb=%f | gc=%f', B, M, gb, gc));
    
      
    