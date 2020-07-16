function strfData = get_all_strfs(cellNames, preprocTypes, stimTypes)

    strfData = cell(length(cellNames), 1);
    
    for j = 1:length(cellNames)

        cellName = cellNames{j};
        cellDir = fullfile('/auto/fdata/mschachter/data', cellName);
        cdata = struct;
        cdata.cellName = cellName;

        for k = 1:length(preprocTypes)

            preprocType = preprocTypes{k};            
    
            for m = 1:length(stimTypes)
                
                stimType = stimTypes{m};
                outputDir = fullfile(cellDir, stimType, 'output');
                
                outputFileLin = sprintf('strflab.tfType_%s.%s.mat', preprocType, stimType);
                outputFile = fullfile(outputDir, outputFileLin);

                vars = load(outputFile);

                strf = vars.modelParamsTrained.w1;
                bias = vars.modelParamsTrained.b1;

                cdata.(preprocType).(stimType).strf = strf;
                cdata.(preprocType).(stimType).bias = bias;
    
            end

        end
        
        strfData{j} = cdata;
        
    end
    
    