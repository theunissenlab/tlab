function run_staggard_nl_fit(cellDir, stimsDir, strfLength, trainingGroups, validationGroups, tfType, tfParams, preprocDesc, stimType, outputFileName)

    %% preprocess stimulus and response
    zscore.stim = 1;
    zscore.resp = 0;
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
        surpriseParams.outputDesc = 'stft.default';
        [surpriseStimLouder, surpriseStimQuieter, groupIndex, surpriseStimInfo, surpriseParams] = preprocSoundSurprise(wholeStim, groupIndex, surpriseParams);
        wholeStim = [surpriseStimLouder; surpriseStimQuieter]';
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

    
    %% get training and early stopping dataset indicies
    trainingIndex = findIdx(trainingGroups, groupIndex);
    validationIndex = findIdx(validationGroups, groupIndex);

    %% Initialize strflab global variables with our stim and responses
    global globDat
    strfData(wholeStim, wholeResponse, groupIndex);

    %% Initialize a linear model
    strfDelays = 0:(strfLength-1);
    linModel = linInit(stimInfo.numStimFeatures, strfDelays, 'linear');
    linModel.b1 = mean(wholeResponse);
    
    %% Initialize a nonlinear model
    %nlModel = glogInit();
    nlModel = gompInit(1);
    
    %% Initialize a linear-nonlinear model
    lnlModel = lnlInit(linModel, nlModel);
    
    %% Initialize linear optimization options
    linOpts = trnDirectFit();
    linOpts.tolerances = [1e-03 5e-04 1e-04 5e-05];
    linOpts.sparsenesses = [0 1 2 4 6 8 10 12];
    linOpts.infoFreqCutoff = 90;
    linOpts.infoWindowSize = 0.250;
    
    %% Initialize nonlinear optimization routines
    nlOpts = trnSCG();
    
    %% Initialze staggard optimization routines
    optOptions = trnStaggardNL();
    optOptions.jackknife = 1;
    optOptions.maxIter = 25;
    optOptions.linOptOptions = linOpts;
    optOptions.nlOptOptions = nlOpts;
    optOptions.sampleRate = stimInfo.sampleRate;
    optOptions.infoWindowSize = 0.250;
    optOptions.infoFreqCutoff = 90;
    
    %% Run staggard optimization
    [modelParamsTrained, optOptions] = strfOpt(lnlModel, trainingIndex, optOptions, validationIndex);
    
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
    