function run_strflab_directfit(cellDir, stimsDir, strfLength, trainingGroups, validationGroups, tfType, tfParams, preprocDesc, stimType, outputFileName)

    tic;

    %% preprocess stimulus and response
    zscore = 1;
    preprocData = preproc_stim_response(cellDir, stimsDir, zscore, stimType, tfType, tfParams, preprocDesc);    
    groupIndex = preprocData.groupIndex;
    wholeStim = preprocData.stim;
    wholeResponse = preprocData.resp;
    dataDir = preprocData.dataDir;
    stimInfo = preprocData.stimInfo;
    
    %% compute surprise if requested
    if isfield(tfParams, 'surprise') && tfParams.surprise        
        fprintf('Computing surprise of spectrogram...\n');
        surpriseParams = struct;
        surpriseParams.outputPath = fullfile(cellDir, stimType, 'preproc');
        surpriseParams.outputDesc = 'default';
        [surpriseStimLouder, surpriseStimQuieter, groupIndex, surpriseStimInfo, surpriseParams] = preprocSoundSurprise(wholeStim, groupIndex, surpriseParams);
        wholeStim = [surpriseStimLouder; surpriseStimQuieter]';
        newStimSize = size(wholeStim)
    end
    
    
    %% check output directory to see if this has already been run
    outputDir = fullfile(dataDir, 'output');
    [s,mess,messid] = mkdir(outputDir);  
    outputFilePath = fullfile(outputDir, outputFileName);

    overwrite = 1;
    if ~overwrite
        fid = fopen(outputFilePath);
        if fid ~= -1
            fclose(fid);
            fprintf('Output file already exists for parameters, returning...\n');
            return;
        end
    end
    
    displayStimResp = 0;
    if displayStimResp       
        display_stim_resp(preprocData);
    end
    
    %% get training and early stopping dataset indicies
    trainingIndex = findIdx(trainingGroups, groupIndex);
    validationIndex = findIdx(validationGroups, groupIndex);
    
    
    %% Initialize strflab global variables with our stim and responses
    global globDat
    strfData(wholeStim, wholeResponse, groupIndex);


    %% Initialize a linear model
    strfDelays = 0:(strfLength-1);
    modelParams = linInit(stimInfo.numStimFeatures, strfDelays, 'linear');
    modelParams.b1 = 0;
    
    
    %% Initialize Direct Fit options
    optOptions = trnDirectFit();
    optOptions.tolerances = [0.100 0.050 0.010 0.005 1e-03 5e-04 1e-04 5e-05];
    optOptions.sparsenesses = [0 1 2 4 8 16 24];
    optOptions.infoFreqCutoff = 90;
    optOptions.infoWindowSize = 0.250;
       
    %% run direct fit
    [modelParamsTrained, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);
    modelParams.b1 = 0;
       
    perfData = compute_performance(preprocData, trainingIndex, validationIndex, modelParamsTrained);
    
    elapsedTime = toc;
    
    fprintf('\nTook %0.0f minutes.\n', elapsedTime/60);
    
    %% save stuff
    save(outputFilePath, 'modelParamsTrained', 'preprocData', 'perfData', ...
                         'trainingIndex', 'validationIndex', ...
                         'cellDir', 'stimsDir', 'strfLength', 'trainingGroups', 'validationGroups', ...
                         'optOptions', 'elapsedTime');
                    
    fprintf('Wrote output file to %s\n', outputFilePath);
                     
    outputFilePathText = strrep(outputFilePath, '.mat', '.txt');
    fid = fopen(outputFilePathText, 'w');
    fprintf(fid, '%d,%f,%f\n', -1, perfData.train.perf, perfData.valid.perf);
    fclose(fid);
                     
    %% display stuff
    display = 1;
    if display
        display_perfdata(preprocData, perfData, modelParamsTrained, trainingIndex, validationIndex);
    end
     %}