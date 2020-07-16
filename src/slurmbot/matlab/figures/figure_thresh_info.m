function figure_thresh_info(cellName)

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
    
    
    types = {'df', 'tg', 'lars'};
    
    thresholds = 0:0.05:0.3;
    
    odata = struct;
    
    for k = 1:length(types)

        t = types{k};
        mResp = mdata.(t).linear.respValid;
        mResp(isnan(mResp)) = 0;
        
        t = types{k};
        odata.(t) = struct;
        odata.(t).resps = zeros(length(thresholds), length(mResp));
        odata.(t).infos = zeros(length(thresholds), 1);
        odata.(t).perfs = zeros(length(thresholds), 1);
        
        for m = 1:length(thresholds)            
            tmResp = mResp;

            tmResp(tmResp < thresholds(m)) = 0;
            [cBoundValid, cValid] = compute_coherence_full(tmResp, respReal, respRealH1, respRealH2, mdata.stimInfo.sampleRate, mdata.avgNumTrials, infoFreqCutoff, infoWindowSize);
            
            odata.(t).resps(m, :) = tmResp;
            odata.(t).infos(m) = cValid.info;
            odata.(t).perfs(m) = cValid.info / cBoundValid.info;            
        end
        
    end
    
    
    for k = 1:length(types)
        
        t = types{k};
        figure('Name', t); hold on;
        for m = 1:length(thresholds)
           
            tmResp = odata.(t).resps(m, :);
            info = odata.(t).infos(m);
            perf = odata.(t).perfs(m);
            
            subplot(length(thresholds), 1, m); hold on;
            plot(respReal, 'k-');
            plot(tmResp, 'r-');
            t1 = sprintf('thresh=%0.2f | info=%f | perf=%f', thresholds(m), info, perf);
            title(t1);
            
        end        
    end
    