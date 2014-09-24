function run_strflab_threshgrad(cellDir, stimsDir, strfLength, nlType, trainingGroups, validationGroups, hyperParams, fileDesc, maxIters, earlyStop)

    if nargin < 9
        maxIters = 8500;
    end
    
    if nargin < 10
        earlyStop = 1;
    end

    tic;
  
    preprocData = preproc_stim_response(cellDir, stimsDir);    
    groupIndex = preprocData.groupIndex;
    wholeStim = preprocData.stim;
    wholeResponse = preprocData.resp;
    dataDir = preprocData.dataDir;
    stimInfo = preprocData.stimInfo;
    
    %% check output directory to see if this has already been run
    outputDir = fullfile(dataDir, 'output');
    [s,mess,messid] = mkdir(outputDir);  
    outputFileName = sprintf('strflab.stft.directfit.%s.mat', fileDesc);
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
    modelParams = linInit(stimInfo.numStimFeatures, strfDelays, nlType);
    if strcmp('linear', nlType)
        modelParams.b1 = mean(wholeResponse);
    elseif strcmp('exponential', nlType)
        wr = wholeResponse;
        wr(wr <= 0) = 1e-3;
        modelParams.b1 = mean(log(wr));
    end
    
    %% Initialize Threshold GradDesc options    
    optOptions = trnThreshGradDesc();
    optOptions.threshold = hyperParams.threshold;
    optOptions.maxIter = maxIters;    
    optOptions.display = 1;
    optOptions.earlyStop = earlyStop;    
    
    %% find optimal step size
    th = optOptions.threshold;
    if strcmp('exponential', nlType)
        if th >= 0 && th < 0.25
            stepSize = 1e-3;
        elseif th >= 0.25 && th < 0.50
            stepSize = 3e-3;
        elseif th >= 0.50 && th < 0.75
            stepSize = 5e-3;
        elseif th >= 0.75 && th < 1.0
            stepSize = 1e-2;
        else
            stepSize = 1;
        end
    else
        if th >= 0 && th < 0.25
            stepSize = 2e-5;
        elseif th >= 0.25 && th < 0.50
            stepSize = 3e-5;
        elseif th >= 0.50 && th < 0.75
            stepSize = 4e-5;
        elseif th >= 0.75 && th < 1.0
            stepSize = 1e-4;
        else
            stepSize = 1e-2;
        end
    end
    stepSize
    
    optOptions.stepSize = stepSize;
    
    %% run threshold gradient descent
    if earlyStop
        [modelParamsTrained, optOptions] = strfOpt(modelParams, trainingIndex, optOptions, validationIndex);
    else
        [modelParamsTrained, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);
    end
   
    
    perfData = compute_performance(preprocData, trainingIndex, validationIndex, modelParamsTrained);
    
    elapsedTime = toc;
    
    fprintf('\nTook %0.0f minutes.\n', elapsedTime/60);
    
    %% save stuff
    save(outputFilePath, 'modelParamsTrained', 'preprocData', 'perfData', ...
                         'trainingIndex', 'validationIndex', ...
                         'cellDir', 'stimsDir', 'strfLength', 'trainingGroups', 'validationGroups', 'fileDesc', ...
                         'optOptions', 'elapsedTime');
                    
    outputFilePathText = strrep(outputFilePath, '.mat', '.txt');
    fid = fopen(outputFilePathText, 'w');
    fprintf(fid, '%d,%f,%f\n', -1, perfData.train.perf, perfData.valid.perf);
    fclose(fid);
                     
    %% display stuff
    display = 1;
    if display
        display_perfdata(preprocData, perfData, modelParamsTrained, trainingIndex, validationIndex);
    end
        