function batch_coherence_visualize(rootdir, stimdir)

    outdir = '~/text/docs/berkeley/coherence_datasets';

    %% collect data
    datasets = find_datasets(rootdir, stimdir, 'conspecific|zfsongs');
    
    infoVals = zeros(length(datasets), 1);
    spikeRate = zeros(length(datasets), 1);    
    cvals = cell(length(datasets), 1);

    for k = 1:length(datasets)
        
        ds = datasets{k};        
        dataDir = ds.dirname;
        
        ifile = fullfile(dataDir, 'info', 'cbnd.mat');
        ivars = load(ifile);
        cbnd = ivars.cbnd;
        
        cvals{k} = cbnd;
        infoVals(k) = cbnd.info;        
        spikeRate(k) = ivars.spikeRate;
    end
    
    goodIndx = find(~isnan(infoVals));
    infoVals = infoVals(goodIndx);
    spikeRate = spikeRate(goodIndx);
    cvals = cvals(goodIndx);
    datasets = datasets(goodIndx);

    infoMean = mean(infoVals);
    infoMode = mode(infoVals);
    spikeMean = mean(spikeRate);
    spikeMode = mode(spikeRate);
    
    [sivals, sindx] = sort(infoVals);
    maxInfo = max(infoVals);
    
    %% plot histograms
    f = figure; hold on;
    hist(infoVals, 25);
    xlabel('NMI (bits)');
    ylabel('# of datasets');
    title(sprintf('Distribution of Info Values: mean=%f, mode=%f', infoMean, infoMode));
    fname = fullfile(outdir, 'nmi_hist.svg');
    plot2svg_2d(fname, f);
        
    f = figure; hold on;
    hist(spikeRate, 25);
    xlabel('Rate (spikes/s)');
    ylabel('# of datasets');
    title(sprintf('Distribution of Spike Rates: mean=%f, mode=%f', spikeMean, spikeMode));
    fname = fullfile(outdir, 'spike_hist.svg');
    plot2svg_2d(fname, f);
    
    %% plot top 10
    topTenIndx = sindx((end-10):end);
    lstrs = cell(length(topTenIndx), 1);
    f = figure; hold on;
    for k = 1:length(topTenIndx)
        
        kindx = topTenIndx(k);
        cbnd = cvals{kindx};
        ds = datasets{kindx};
        
        [pathstr, cellName, ext, versn] = fileparts(fileparts(ds.dirname));
        [pathstr, groupName, ext, versn] = fileparts(pathstr);
        lstrs{k} = [groupName '.' cellName]; 
        
        ninfo = cbnd.info / maxInfo;
        clr = [ninfo 0 (1-ninfo)];
        plot(cbnd.f, cbnd.c, '-', 'Color', clr, 'LineWidth', 2);
    end
    xlabel('Frequency (Hz)');
    ylabel('Coherence');
    title('Top 10 Datasets');
    legend(char(lstrs));    
    fname = fullfile(outdir, 'top10_coherence.svg');
    plot2svg_2d(fname, f);
    
    %% plot worst 10
    bottomTenIndx = sindx(1:10);
    lstrs = cell(length(bottomTenIndx), 1);
    f = figure; hold on;
    for k = 1:length(bottomTenIndx)
        
        kindx = bottomTenIndx(k);
        cbnd = cvals{kindx};
        ds = datasets{kindx};
        
        [pathstr, cellName, ext, versn] = fileparts(fileparts(ds.dirname));
        [pathstr, groupName, ext, versn] = fileparts(pathstr);
        lstrs{k} = [groupName '.' cellName]; 
        
        ninfo = cbnd.info / maxInfo;
        clr = [ninfo 0 (1-ninfo)];
        plot(cbnd.f, cbnd.c, '-', 'Color', clr, 'LineWidth', 2);
    end
    title('Worst 10 Datasets');
    legend(char(lstrs));    
    fname = fullfile(outdir, 'worst10_coherence.svg');
    plot2svg_2d(fname, f);

    %% write all to file
    fpath = fullfile(outdir, 'conspecific+zfsongs_coherence.csv');
    fid = fopen(fpath, 'w');
    for k = 1:length(sindx)
        kindx = sindx(k);
        ds = datasets{kindx};
        
        [pathstr, cellName, ext, versn] = fileparts(fileparts(ds.dirname));
        [pathstr, groupName, ext, versn] = fileparts(pathstr);
        dsName = [groupName '.' cellName]; 
        
        fprintf(fid, '%s,%f,%f\n', dsName, infoVals(kindx), spikeRate(kindx));                
    end
    fclose(fid);
