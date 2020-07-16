function fig_compare_df_threshgrad(mdatas)

    nParts = size(mdatas, 1);
    
    dfPerfs = zeros(nParts, 1);
    tgPerfs = zeros(nParts, 1);
    
    figure; hold on;
    for k = 1:nParts
        
        dfData = mdatas{k, 1};
        tgData = mdatas{k, 2};
        
        dfStrf = dfData.model.strf;
        tgStrf = tgData.model.strf;
        
        dfInfo = dfData.coherence.validation.info_mean;
        tgInfo = tgData.coherence.validation.info_mean;
        
        dfPerfs(k) = dfInfo;
        tgPerfs(k) = tgInfo;
        
        dfIndx = (k-1)*2 + 1;
        tgIndx = dfIndx + 1;
        
        subplot(nParts, 2, dfIndx); hold on;
        imagesc(dfStrf); axis tight;
        cval = max(abs(dfStrf(:)));
        caxis([-cval cval]);
        title(sprintf('DF: info=%0.1f bits', dfInfo));
        
        subplot(nParts, 2, tgIndx); hold on;
        imagesc(tgStrf); axis tight;
        cval = max(abs(tgStrf(:)));
        caxis([-cval cval]);
        title(sprintf('TG: info=%0.1f bits', tgInfo));
    end
    
    fprintf('DF Perf: %0.1f +/- %0.2f\n', mean(dfPerfs), std(dfPerfs));
    fprintf('TG Perf: %0.1f +/- %0.2f\n', mean(tgPerfs), std(tgPerfs));
    
    pdiff = tgPerfs - dfPerfs;    
    fprintf('Diff in Perf: %0.1f +/- %0.2f\n', mean(pdiff), std(pdiff));
    