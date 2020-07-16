function figure_response_hist(cellName)

    dataDir = '/auto/fdata/mschachter/data';
    stimsDir = fullfile(dataDir, 'all_stims');

    mdata = get_method_data(dataDir, stimsDir, cellName);
    
    mTypes = {'df', 'tg', 'lars'};
    nlTypes = {'linear', 'exponential'};
    
    nhist = 50;
    cutoff = 3e-1;
    
    histFig = figure('Name', cellName); hold on;
    
    subplot(4, 2, 1);
    resp = mdata.resp;
    resp(resp <= cutoff) = [];
    hist(resp, nhist);
    title('Real');
            
    vindx = findIdx(mdata.vGroups, mdata.groupIndex);
    respRealValid = mdata.resp(vindx);
    
    subplot(4, 2, 2);
    hist(resp, nhist);
    title('Real');    
    
    for m = 1:length(nlTypes)
        
        for k = 1:length(mTypes)
        
            tdata = mdata.(mTypes{k}).(nlTypes{m});
            
            %% raw response data
            figure; hold on;
            plot(respRealValid, 'k-', 'LineWidth', 2);
            plot(tdata.respValid, 'r-');
            tname = sprintf('%s | %s | %0.4f', nlTypes{m}, mTypes{k}, tdata.perf);
            title(tname); 
            
            %% histogram data
            r = [tdata.respTrain' tdata.respValid'];
            r(r <= cutoff) = [];
            figure(histFig);
            sp = 2 + 2*(k-1) + 1 + mod(m+1, 2);
            subplot(4, 2, sp);
            hist(r, nhist);
            
            tname = sprintf('%s | %s | %0.4f', nlTypes{m}, mTypes{k}, tdata.perf);
            title(tname);
        end
    end    
    