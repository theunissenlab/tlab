function fit_isi_dataset(rootDir, stimDir, cellDirs)

    dirNames = cell(length(cellDirs), 1);
    for k = 1:length(cellDirs)
        dirNames{k} = fullfile(rootDir, cellDirs{k});
    end
        
    preprocType = 'ft';
    
    for k = 1:length(dirNames)
        
        datasets = find_datasets(dirNames{k}, stimDir, 'conspecific|zfsongs');
        srFiles = datasets{1}.srPairs;
        dataDir = datasets{1}.dirname;

        %% preprocess stimuli and responses
        stimFiles = srFiles.stimFiles;
        respFiles = srFiles.respFiles;

        preprocDir = fullfile(dataDir, 'preproc');
        [s,mess,messid] = mkdir(preprocDir);  

        stimOutPattern = ['stim.preproc-' preprocType '.%d.mat'];
        respOutPattern = ['resp.preproc-' preprocType '.%d.mat'];

        srData = preprocess_sound(stimFiles, respFiles, preprocType, struct, preprocDir, stimOutPattern, respOutPattern);
        
        fit_isi(srData, cellDirs{k});
        
    end