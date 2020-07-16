function display_strf_data(outputFileName)

    vars = load(outputFileName);
    
    modelParamsTrained = vars.modelParamsTrained;
    preprocData = vars.preprocData;
    perfData = vars.perfData;
    trainingIndex = vars.trainingIndex;
    validationIndex = vars.validationIndex;
    
    trainingGroups = vars.trainingGroups
    validationGroups = vars.validationGroups
    
    clear vars;
    
    display_stim_resp(preprocData);    
    display_perfdata(preprocData, perfData, modelParamsTrained, trainingIndex, validationIndex);
    