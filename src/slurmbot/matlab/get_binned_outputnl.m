function [xvals, yvals, stds] = get_binned_outputnl(linResp, psth, convWithGauss)

    if nargin < 3
        convWithGauss = 0;
    end
    
    if convWithGauss        
        gwidth = 0.003;
        grng = -0.05:0.001:0.05;
        gconv = exp(-(grng / gwidth).^2);
        
        spsth = conv(psth, gconv, 'same');
        spsth = (spsth / max(spsth))*max(psth);
        
        psth = spsth;
    end
    

    [nLinVals, linBinVals] = hist(linResp, 35);
    psthDists = cell(length(linBinVals)-1, 1);
    psthMeans = zeros(length(linBinVals)-1, 1);
    psthStds = zeros(length(linBinVals)-1, 1);
    for k = 2:length(linBinVals)        
        lowval = linBinVals(k-1);
        highval = linBinVals(k);
        pdist = psth((linResp > lowval) & (linResp <= highval));
        psthMeans(k-1) = mean(pdist);
        psthStds(k-1) = std(pdist);
        psthDists{k-1} = pdist;        
    end
    
    xvals = cv(linBinVals(2:end));
    yvals = cv(psthMeans);
    stds = cv(psthStds);
    