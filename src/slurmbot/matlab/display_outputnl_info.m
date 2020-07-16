function [xvals, yvals, stds, perfData, gompParams] = display_outputnl_info(outputFileName, svgFile)

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
    
    %psthThresh = 0.2;
    %respTrain(respTrain < psthThresh) = 0;
    
    %% get psth-conditioned distributions    
    psthVals = unique(realResp);
    linDists = cell(length(psthVals), 1);
    linMeans = zeros(length(psthVals), 1);
    linStds = zeros(length(psthVals), 1);
    for k = 1:length(psthVals)               
        ldist = respTrain(realResp == psthVals(k));        
        lmean = mean(ldist);
        linMeans(k) = lmean;
        linDists{k} = ldist;    
        linStds(k) = std(ldist);
    end
    
    %{
    figure('Name', 'P(X)'); hold on;
    hist(respTrain, 75);
    
    figure('Name', 'P(Y)'); hold on;
    hist(realResp, 20);
    %}
    %{
    figure('Name', 'P(X|Y)');hold on;
    for k = 1:length(psthVals)        
       ldist = linDists{k};       
       subplot(length(psthVals), 1, k); hold on;
       hist(ldist, 25);
       %title(sprintf('psth=%0.2f | mean=%0.3f', psthVals(k), lmean));
    end
    
    figure('Name', 'E[P(X|Y)]');    
    errorbar(linMeans, psthVals, linStds);
    %}
    
    %% get linear-response-conditioned distributions
    
    [nLinVals, linBinVals] = hist(respTrain, 10);
    psthDists = cell(length(linBinVals)-1, 1);
    psthMeans = zeros(length(linBinVals)-1, 1);
    psthStds = zeros(length(linBinVals)-1, 1);
    for k = 2:length(linBinVals)        
        lowval = linBinVals(k-1);
        highval = linBinVals(k);
        pdist = realResp((respTrain > lowval) & (respTrain <= highval));
        psthMeans(k-1) = mean(pdist);
        psthStds(k-1) = std(pdist);
        psthDists{k-1} = pdist;        
    end
    %{
    figure('Name', 'P(Y|X)'); hold on;
    for k = 1:length(psthDists)       
        pdist = psthDists{k};
        subplot(length(psthDists), 1, k);
        hist(pdist, 10);        
    end
    %}
    
    xvals = linBinVals(2:end);
    yvals = psthMeans;
    stds = psthStds;
    
    a = max(yvals);
    [b, c] = fit_gompertz_fmin(cv(xvals), cv(yvals), a, -abs(randn()), -abs(randn()));
    gompParams = [a, b, c];
    
    xinc = (max(xvals) - min(xvals)) / 100;
    gx = min(xvals):xinc:max(xvals);
    
    gy = gompertz(a, b, c, gx);
    h1 = figure('Name', 'E[P(Y|X)]'); hold on;
    errorbar(linBinVals(2:end), psthMeans, psthStds, 'k-', 'LineWidth', 2);
    plot(gx, gy, 'r-', 'LineWidth', 2);
    xlabel('x');
    ylabel('E[P(Y|X=x)]');
    legend('E[P(Y|X)]', 'Gompertz');
    axis tight;
    
    if nargin > 1
        plot2svg_2d(svgFile, h1);
    end
    
    
    
    
    
    
    