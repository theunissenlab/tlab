function modelDatas = get_df_tg_modeldatas(cellName, threshold, preprocType)

    if nargin < 3
        preprocType = 'stft';
    end

    stimClass = 'Con';
    
    switch preprocType        
        case 'stft'
            preprocDesc = 'stft.nstd_6.fband_125';
        case 'lyons'
            [earQ, step] = get_best_lyons_params(cellName, stimClass);
            preprocDesc = sprintf('lyons.agc_1.earQ_%d.step_%0.2f', earQ, step);
        case 'surprise'
            preprocDesc = sprintf('surprise.dfw_3.dtw_3.dg_4.%s', stimClass);
    end
    
    dfDesc = sprintf('directfit.%s.%s.%%d.h5', stimClass, preprocDesc);
    tgDesc = sprintf('threshgrad.%s.%s.%0.2f.%%d.h5', stimClass, preprocDesc, threshold);
    
    outputDir = fullfile('/auto/k6/mschachter/pystrfs/units', cellName, 'output');
    
    numParts = 0:9;
    
    modelDatas = cell(length(numParts), 2);
    for k = 1:length(numParts)
                
        outputName = fullfile(outputDir, sprintf(dfDesc, numParts(k)));
        mdata_df = get_model_data(outputName);
        
        outputName = fullfile(outputDir, sprintf(tgDesc, numParts(k)));
        mdata_tg = get_model_data(outputName);
        
        modelDatas{k, 1} = mdata_df;
        modelDatas{k, 2} = mdata_tg;
        
    end
    