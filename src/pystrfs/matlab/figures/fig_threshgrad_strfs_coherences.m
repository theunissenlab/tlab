function fig_threshgrad_strfs_coherences(unit)

    outputDir = sprintf('/auto/k6/mschachter/pystrfs/units/%s/output', unit);
    fpat = 'threshgrad.linear.Con.stft.nstd_6.fband_125.thresh_0.75.0.h5';
    
    outputFile = fullfile(outputDir, fpat); 
    
    mdata = get_model_data(outputFile);
    
    cb = mdata.coherence.bound;
    cm = mdata.coherence.validation;
    
    figure; hold on;
    plot(cb.f, cb.upper, 'r-', 'LineWidth', 3);
    plot(cb.f, cb.mean, 'k-', 'LineWidth', 3);
    plot(cb.f, cb.lower, 'b-', 'LineWidth', 3);
    
    plot(cm.f, cm.upper, 'r--', 'LineWidth', 3);
    plot(cm.f, cm.mean, 'k--', 'LineWidth', 3);
    plot(cm.f, cm.lower, 'b--', 'LineWidth', 3);
    
    axis tight;
    