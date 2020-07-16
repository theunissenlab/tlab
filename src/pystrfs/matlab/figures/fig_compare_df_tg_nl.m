function fig_compare_df_tg_nl(cellName, preprocType, stimClass)

    data = get_df_tg_nl_cube(cellName, preprocType, stimClass);    
    infos = data.infos;
    strfs = data.strfs;
    
    figName = sprintf('%s | %s | %s', cellName, preprocType, stimClass);
    fprintf('\n%s\n', figName);
    fprintf('--------------------------\n');
    
    %% compare df to threshgrad
    dfPerf = infos(1, :, 2) ./ infos(1, :, 1);
    tgPerf = infos(4, :, 2) ./ infos(4, :, 1);
    df_tg_diffs = tgPerf - dfPerf;
    
    fprintf('  DF vs. Threshgrad\n');
    fprintf('    DF Performance: %0.3f +/- %0.3f\n', mean(dfPerf), std(dfPerf));
    fprintf('    TG Performance: %0.3f +/- %0.3f\n', mean(tgPerf), std(tgPerf));
    fprintf('    Difference: %0.4f +/- %0.4f\n\n', mean(df_tg_diffs), std(df_tg_diffs));
    
    %% compare df+nl (dist) to tg
    dfnlPerf = infos(2, :, 2) ./ infos(2, :, 1);
    tgPerf = infos(4, :, 2) ./ infos(4, :, 1);
    dfnl_tg_diffs = tgPerf - dfnlPerf;
    
    fprintf('  DF+NL (dist) vs. Threshgrad\n');
    fprintf('    DF+NL(dist) Performance: %0.3f +/- %0.3f\n', mean(dfnlPerf), std(dfnlPerf));
    fprintf('    TG+exp Performance: %0.3f +/- %0.3f\n', mean(tgPerf), std(tgPerf));
    fprintf('    Difference: %0.4f +/- %0.4f\n\n', mean(dfnl_tg_diffs), std(dfnl_tg_diffs));
    
    %% compare df+nl (spline) to tg
    dfnlPerf = infos(3, :, 2) ./ infos(3, :, 1);
    tgPerf = infos(4, :, 2) ./ infos(4, :, 1);
    dfnl_tg_diffs = tgPerf - dfnlPerf;
    
    fprintf('  DF+NL (spline) vs. Threshgrad\n');
    fprintf('    DF+NL(spline) Performance: %0.3f +/- %0.3f\n', mean(dfnlPerf), std(dfnlPerf));
    fprintf('    TG+exp Performance: %0.3f +/- %0.3f\n', mean(tgPerf), std(tgPerf));
    fprintf('    Difference: %0.4f +/- %0.4f\n\n', mean(dfnl_tg_diffs), std(dfnl_tg_diffs));
    
    %% compare df+nl to exp (dist)
    dfnlPerf = infos(2, :, 2) ./ infos(2, :, 1);
    tgexpPerf = infos(5, :, 2) ./ infos(5, :, 1);
    dfnl_tgexp_diffs = tgexpPerf - dfnlPerf;
    
    fprintf('  DF+NL (dist) vs. Threshgrad + exp\n');
    fprintf('    DF+NL(dist) Performance: %0.3f +/- %0.3f\n', mean(dfnlPerf), std(dfnlPerf));
    fprintf('    TG+exp Performance: %0.3f +/- %0.3f\n', mean(tgexpPerf), std(tgexpPerf));
    fprintf('    Difference: %0.4f +/- %0.4f\n\n', mean(dfnl_tgexp_diffs), std(dfnl_tgexp_diffs));
    
    %% compare df+nl to exp (spline)
    dfnlPerf = infos(3, :, 2) ./ infos(3, :, 1);
    tgexpPerf = infos(5, :, 2) ./ infos(5, :, 1);
    dfnl_tgexp_diffs = tgexpPerf - dfnlPerf;
    
    fprintf('  DF+NL (spline) vs. Threshgrad + exp\n');
    fprintf('    DF+NL(spline) Performance: %0.3f +/- %0.3f\n', mean(dfnlPerf), std(dfnlPerf));
    fprintf('    TG+exp Performance: %0.3f +/- %0.3f\n', mean(tgexpPerf), std(tgexpPerf));
    fprintf('    Difference: %0.4f +/- %0.4f\n\n', mean(dfnl_tgexp_diffs), std(dfnl_tgexp_diffs));
    
    %% compare df+nl to log (dist)
    dfnlPerf = infos(2, :, 2) ./ infos(2, :, 1);
    tglogPerf = infos(6, :, 2) ./ infos(6, :, 1);
    dfnl_tglog_diffs = tglogPerf - dfnlPerf;
    
    fprintf('  DF+NL (dist) vs. Threshgrad + log\n');
    fprintf('    DF+NL(dist) Performance: %0.3f +/- %0.3f\n', mean(dfnlPerf), std(dfnlPerf));
    fprintf('    TG+log Performance: %0.3f +/- %0.3f\n', mean(tglogPerf), std(tglogPerf));
    fprintf('    Difference: %0.4f +/- %0.4f\n\n', mean(dfnl_tglog_diffs), std(dfnl_tglog_diffs));
    
    %% compare df+nl to log (spline)
    dfnlPerf = infos(3, :, 2) ./ infos(3, :, 1);
    tglogPerf = infos(6, :, 2) ./ infos(6, :, 1);
    dfnl_tglog_diffs = tglogPerf - dfnlPerf;
    
    fprintf('  DF+NL (spline) vs. Threshgrad + log\n');
    fprintf('    DF+NL(spline) Performance: %0.3f +/- %0.3f\n', mean(dfnlPerf), std(dfnlPerf));
    fprintf('    TG+log Performance: %0.3f +/- %0.3f\n', mean(tglogPerf), std(tglogPerf));
    fprintf('    Difference: %0.4f +/- %0.4f\n\n', mean(dfnl_tglog_diffs), std(dfnl_tglog_diffs));
    
    %% compare exp to log    
    tgexpPerf = infos(5, :, 2) ./ infos(5, :, 1);
    tglogPerf = infos(6, :, 2) ./ infos(6, :, 1);
    tgexp_tglog_diffs = tglogPerf - tgexpPerf;
    
    fprintf('  Threshgrad+exp vs. Threshgrad+log\n');
    fprintf('    TG+exp Performance: %0.3f +/- %0.3f\n', mean(tgexpPerf), std(tgexpPerf));
    fprintf('    TG+log Performance: %0.3f +/- %0.3f\n', mean(tglogPerf), std(tglogPerf));
    fprintf('    Difference: %0.4f +/- %0.4f\n\n', mean(tgexp_tglog_diffs), std(tgexp_tglog_diffs));
    
    
    %% visualize strfs of TG, exp, log
    nprts = size(strfs, 2);
    figure('Name', figName); hold on;
    for k = 1:nprts       
        indx = (k-1)*3 + 1;        
        tgStrf = strfs{4, k};
        tgexpStrf = strfs{5, k};
        tglogStrf = strfs{6, k};
        
        subplot(nprts, 3, indx); hold on;
        cval = max(abs(tgStrf(:)));
        imagesc(tgStrf); axis tight;
        caxis([-cval cval]);
        
        subplot(nprts, 3, indx+1); hold on;
        cval = max(abs(tgexpStrf(:)));
        imagesc(tgexpStrf); axis tight;
        caxis([-cval cval]);
        
        subplot(nprts, 3, indx+2); hold on;
        cval = max(abs(tglogStrf(:)));
        imagesc(tglogStrf); axis tight;
        caxis([-cval cval]);                
    end
    
    
    
    
    
    
    
    
    
    
    
    