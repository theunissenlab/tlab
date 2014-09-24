function fig_threshgrad_strfs_thresholds(unit)

    outputDir = sprintf('/auto/k6/mschachter/pystrfs/units/%s/output', unit);
    fpat = 'threshgrad.linear.Con.stft.nstd_6.fband_125.thresh_%0.2f.0.h5';
    
    thresholds = [0.0, 0.25, 0.5, 0.75, 1.0];
        
    for k = 1:length(thresholds)        
        thresh = thresholds(k);    
        outputFile = sprintf(fullfile(outputDir, fpat), thresh);       
        mdata = get_model_data(outputFile);
        
        figure; hold on;
        imagesc(mdata.model.strf); axis tight;
        title(sprintf('Threshold=%0.2f', thresh));
    end
    