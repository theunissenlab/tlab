function figure_thresh_info2(cellName)

    dataRootDir = '/auto/fdata/mschachter/data';
    stimsDir = fullfile(dataRootDir, 'all_stims');
    
    infoFreqCutoff = 90; %Hz
    infoWindowSize = 0.250; %250ms
    
    mdata = get_method_data(dataRootDir, stimsDir, cellName);
    
    wholeResp = mdata.resp;
    wholeRespH1 = mdata.respHalf1;
    wholeRespH2 = mdata.respHalf2;
    
    vindx = findIdx(mdata.vGroups, mdata.groupIndex);
    
    respReal = wholeResp(vindx);
    respRealH1 = wholeRespH1(vindx);
    respRealH2 = wholeRespH2(vindx);
    
    thresh = 0.2;    
    respValid = mdata.tg.linear.respValid;
    
    respValidThresh = respValid;
    respValidThresh(respValidThresh < thresh) = 0;
    
    [cBoundValid, cValid] = compute_coherence_full(respValid, respReal, respRealH1, respRealH2, mdata.stimInfo.sampleRate, mdata.avgNumTrials, infoFreqCutoff, infoWindowSize);
            
    [cBoundValidThresh, cValidThresh] = compute_coherence_full(respValidThresh, respReal, respRealH1, respRealH2, mdata.stimInfo.sampleRate, mdata.avgNumTrials, infoFreqCutoff, infoWindowSize);

    perf = cValid.info / cBoundValid.info;
    tperf = cValidThresh.info / cBoundValidThresh.info;    
    
    t = 0:0.001:(length(respReal)-1)/1000;
    h1 = figure; hold on;
    plot(t, respReal, 'k-', 'LineWidth', 2);
    plot(t, respValid, 'r-', 'LineWidth', 2);
    legend('Real PSTH', 'TG Pred.');
    axis tight;
    title(sprintf('perf=%f', perf));
    
    h2 = figure; hold on;
    plot(t, respReal, 'k-', 'LineWidth', 2);
    plot(t, respValidThresh, 'r-', 'LineWidth', 2);
    legend('Real PSTH', 'TG Pred. w/ Thresh NL');
    axis tight;
    title(sprintf('perf=%f', tperf));
    
    ofile1 = '~/Desktop/sfn_poster/images/figure_thresh_pre.svg';
    plot2svg(ofile1, h1);
    ofile2 = '~/Desktop/sfn_poster/images/figure_thresh_post.svg';
    plot2svg(ofile2, h2);
    
    