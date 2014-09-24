function test_threshgrad_outputnl(cellName, preprocName, threshold)

    h5 = h5utils();
    
    %cellName = 'yy1617_4_A';
    cellDir= fullfile('/auto/k6/mschachter/pystrfs/units', cellName);
    unitFile = fullfile(cellDir, 'unit.h5');
    outputDir = fullfile(cellDir, 'output');    
    
    %preprocName = 'stft.nstd_6.fband_125';
    preprocDir = '/auto/k6/mschachter/pystrfs/preproc';
    preprocFile = fullfile(preprocDir, sprintf('%s.h5', preprocName));
    
    stimClass = 'Con';
    
    cvRuns = 0:9;
     
    for k = 1:length(cvRuns)
        
        
        
        
        
        
    end