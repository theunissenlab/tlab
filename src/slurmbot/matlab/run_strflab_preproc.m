function run_strflab_preproc(cellDir, stimsDir, preprocType, preprocDesc, tfParams, stimType)


    %% get datafiles in cell directory
    datasets = find_datasets(cellDir, stimsDir, stimType);    
    stimFiles = datasets{1}.srPairs.stimFiles;
    respFiles = datasets{1}.srPairs.respFiles;
    dataDir = datasets{1}.dirname;
   
    
    %% preprocess the sound to produce spectrograms
    preprocDir = fullfile(dataDir, 'preproc');
    [s,mess,messid] = mkdir(preprocDir);  

    preprocStimParams = struct;      
    preprocStimParams.tfType = preprocType;
    preprocStimParams.tfParams = tfParams;
    preprocStimParams.outputDir = preprocDir;
    preprocStimParams.overwrite = 1;
    preprocStimParams.outputPattern = ['stim.preproc-' preprocDesc '.%d.mat'];
    
    [wholeStim, groupIndex, stimInfo, preprocStimParams] = preprocSound(stimFiles, preprocStimParams);
