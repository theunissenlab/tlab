function figure_strf_diffs(cellName, type)
       
    if nargin < 2
        type = 'tg';
    end

    dataDir = '/auto/fdata/mschachter/data';
    stimsDir = '/auto/fdata/mschachter/data/all_stims';

    mdata = get_method_data(dataDir, stimsDir, cellName);
    
    dfStrf = mdata.df.linear.strf;
    oStrf = mdata.(type).linear.strf;
    
    
    dfResp = mdata.df.linear.respValid;
    oResp = mdata.(type).linear.respValid;
    
    dfl1 = sum(abs(dfStrf(:)));
    dfl2 = sum(abs(dfStrf(:).^2));
    dfPerf = mdata.df.linear.perf;
    
    ol1 = sum(abs(oStrf(:)));
    ol2 = sum(abs(oStrf(:).^2));
    oPerf = mdata.(type).linear.perf;
    
    dfs = dfStrf(:);
    dfsNeg = dfs;
    dfsNeg(dfsNeg > 0) = 0;
    dfsPos = dfs;
    dfsPos(dfsPos < 0) = 0;
    
    dfNegStrf = reshape(dfsNeg, 60, 75);
    dfPosStrf = reshape(dfsPos, 60, 75);
    
    dfPosL1 = sum(abs(dfsPos));
    dfPosL2 = sum(dfsPos.^2);
    dfNegL1 = sum(abs(dfsNeg));
    dfNegL2 = sum(dfsNeg.^2);
    
    os = oStrf(:);
    osNeg = os;
    osNeg(osNeg > 0) = 0;
    osPos = os;
    osPos(osPos < 0) = 0;
    
    oNegStrf = reshape(osNeg, 60, 75);
    oPosStrf = reshape(osPos, 60, 75);    
    
    oPosL1 = sum(abs(osPos));
    oPosL2 = sum(osPos.^2);
    oNegL1 = sum(abs(osNeg));
    oNegL2 = sum(osNeg.^2);
    
    
    vindx = findIdx(mdata.vGroups, mdata.groupIndex);
    realResp = mdata.resp(vindx);
    t = 0:0.001:(length(realResp)-1)/1000;
    h = figure; hold on;    
    plot(t, mdata.df.linear.respValid, 'k-', 'LineWidth', 2);
    plot(t, mdata.(type).linear.respValid, 'r-', 'LineWidth', 2);
    legend('DF', type);
    axis tight;
    xlabel('Time');
    ylabel('P[spike]');
    outputFile = sprintf('~/Desktop/sfn_poster/images/figure_resps_%s.svg', type);
    plot2svg_2d(outputFile, h);
    
    h2 = figure; hold on;
    rdiff = mdata.df.linear.respValid - mdata.(type).linear.respValid;
    plot(t, rdiff, 'k-', 'LineWidth', 2);
    axis tight;
    xlabel('Time');
    ylabel('\DeltaP[spike]');
    outputFile = sprintf('~/Desktop/sfn_poster/images/figure_resp_diff_%s.svg', type);
    plot2svg_2d(outputFile, h2);
    
    
    %% DF vs other
    %{
    figure('Name', sprintf('DF vs %s: %s', type, cellName)); hold on;
    
    subplot(3, 2, 1); hold on;
    smax = max(abs(dfStrf(:)));
    imagesc(dfStrf); axis tight;
    caxis([-smax smax]);
    colorbar;
    t1 = sprintf('l1=%f | l2=%f | perf=%f', dfl1, dfl2, dfPerf);
    title(t1);
    
    subplot(3, 2, 2); hold on;
    plot(dfResp, 'k-'); axis tight;
    
    subplot(3, 2, 3); hold on;
    smax = max(abs(oStrf(:)));
    imagesc(oStrf); axis tight;
    caxis([-smax smax]);
    colorbar;
    t1 = sprintf('l1=%f | l2=%f | perf=%f', ol1, ol2, oPerf);
    title(t1);
    
    subplot(3, 2, 4); hold on;
    plot(oResp, 'k-'); axis tight;
    
    sdiff = dfStrf - oStrf;
    smax = max(abs(sdiff(:)));
    subplot(3, 2, 5); hold on;
    imagesc(sdiff); axis tight;
    caxis([-smax smax]);
    colorbar;
    title('Diff');
    
    rdiff = dfResp - oResp;
    subplot(3, 2, 6); hold on;
    plot(rdiff, 'r-'); axis tight;
    
    
    %% differences in positive and negative parts
    figure('Name', sprintf('Decomp: DF vs %s: %s', type, cellName)); hold on;
    
    subplot(2, 2, 1); hold on;
    imagesc(dfPosStrf); axis tight;
    smax = max(dfPosStrf(:));
    caxis([-smax smax]);
    colorbar;
    t1 = sprintf('Postive | l1=%f | l2=%f', dfPosL1, dfPosL2);
    title(t1);
    
    subplot(2, 2, 2); hold on;
    imagesc(dfNegStrf); axis tight;
    smax = max(abs(dfNegStrf(:)));
    caxis([-smax smax]);
    colorbar;
    t1 = sprintf('Negative | l1=%f | l2=%f', dfNegL1, dfNegL2);
    title(t1);
    
    subplot(2, 2, 3); hold on;
    imagesc(oPosStrf); axis tight;
    smax = max(oPosStrf(:));
    caxis([-smax smax]);
    colorbar;
    t1 = sprintf('Postive | l1=%f | l2=%f', oPosL1, oPosL2);
    title(t1);
    
    subplot(2, 2, 4); hold on;
    imagesc(oNegStrf); axis tight;
    smax = max(abs(oNegStrf(:)));
    caxis([-smax smax]);
    colorbar;
    t1 = sprintf('Negative | l1=%f | l2=%f', oNegL1, oNegL2);
    title(t1);
%}
    
    
    
    
    