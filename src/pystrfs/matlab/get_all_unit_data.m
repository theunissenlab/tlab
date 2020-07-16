
function allData = get_all_unit_data(stimClasses, preprocTypes, unitName)

    allData = struct;    
    
    allData.unitName = unitName;
    unitDir = fullfile('/auto/k6/mschachter/pystrfs/units', unitName);
    outputDir = fullfile(unitDir, 'output');
       
    for k = 1:length(stimClasses)        
        stimClass = stimClasses{k};        
        allData.(stimClass) = struct;        
        for m = 1:length(preprocTypes)            
            preprocType = preprocTypes{m};
            
            switch preprocType
                case 'stft'
                    preprocDesc = 'stft.nstd_6.fband_125';
                case 'lyons'
                    preprocDesc = 'lyons.best';
                case 'surprise'
                    preprocDesc = sprintf('surprise.dfw_3.dtw_3.dg_4.%s', stimClass);
            end
            
            pdata = struct;            
            responseFileBase = sprintf('directfit.%s.%s.h5', stimClass, preprocDesc);
            responseFileDists = sprintf('directfit.nl_dists.%s.%s.h5', stimClass, preprocDesc);
            responseFileSpline = sprintf('directfit.nl_spline.%s.%s.h5', stimClass, preprocDesc);

            mdata = get_model_data(fullfile(outputDir, responseFileBase));
            mdata_dists = get_model_data(fullfile(outputDir, responseFileDists));
            mdata_spline = get_model_data(fullfile(outputDir, responseFileSpline));
            
            pdata.linear = mdata;
            pdata.dists = mdata_dists;
            pdata.spline = mdata_spline;
            
            allData.(stimClass).(preprocType) = pdata;            
        end
    end
   