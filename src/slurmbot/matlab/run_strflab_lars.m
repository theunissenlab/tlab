function run_strflab_lars(cellDir, stimsDir, strfLength, nlType, trainingGroups, validationGroups, hyperParams, fileDesc, maxIters, findBest)

    tic;

    earlyStop = 0;
    if nargin < 8
        maxIters = 750;
    end
    if nargin < 9
        findBest = 1;
    end
    
    preprocData = preproc_stim_response(cellDir, stimsDir);    
    groupIndex = preprocData.groupIndex;
    wholeStim = preprocData.stim;
    wholeResponse = preprocData.resp;    
    dataDir = preprocData.dataDir;
    stimInfo = preprocData.stimInfo;
    
    
    %% check output directory to see if this has already been run
    outputDir = fullfile(dataDir, 'output');
    [s,mess,messid] = mkdir(outputDir);  
    outputFileName = sprintf('strflab.stft.LARS.%s.mat', fileDesc);
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
   
    
    %% train on the inverted nonlinearity if not linear
    wholeResponseNL = wholeResponse;
    if strcmp('exponential', nlType)
        wr = wholeResponse;
        wr(wr <= 0) = 1e-3;
        wholeResponseNL = log(wr);
    end
    
    displayStimResp = 0;
    if displayStimResp       
        display_stim_resp(preprocData);
    end

    %% get training and validation stopping dataset indicies
    trainingIndex = findIdx(trainingGroups, groupIndex);
    validationIndex = findIdx(validationGroups, groupIndex);
    
    
    %% Initialize strflab global variables with our stim and responses
    global globDat
    strfData(wholeStim, wholeResponseNL, groupIndex);


    %% Initialize a linear model
    strfDelays = 0:(strfLength-1);
    modelParams = linInit(stimInfo.numStimFeatures, strfDelays);
    modelParams.b1 = mean(wholeResponseNL);

    
    %% Initialize LARS options
    optOptions = trnLARS();
    optOptions.optName = 'LARS';
    optOptions.type = 'linear';
    optOptions.method = 'en';
    optOptions.lambda2 = hyperParams.lambda2;
    optOptions.display = 1;
    optOptions.earlyStop = earlyStop;
    %optOptions.maxIter = round(0.50*strfLength*stimInfo.numStimFeatures);
    optOptions.maxIter = maxIters;
    
    if earlyStop
        [modelParamsTrained, optOptions] = strfOpt(modelParams, trainingIndex, optOptions, validationIndex);
    else
        [modelParamsTrained, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);
    end
    
    
    %% reset the response data to it's non-inverted form    
    strfData(wholeStim, wholeResponse, groupIndex);
    
    if findBest
        %% go through each of the possible STRFs (one for each iteration), and see how it does on validation sets
        fprintf('Finding best iteration...\n');        
        betas = optOptions.diagnostics.betas;
        bestScore = -1.0;
        bestStrf = [];
        bestIter = -1;
        for k = 1:length(betas)

            strf = reshape(betas{k}, stimInfo.numStimFeatures, strfLength);
            modelParamsCopy = linInit(stimInfo.numStimFeatures, strfDelays, nlType);
            
            modelParamsCopy.b1 = modelParamsTrained.b1;
            modelParamsCopy.w1 = strf;

            perfData = compute_performance(preprocData, [], validationIndex, modelParamsCopy);
            
            if perfData.valid.perf > bestScore
                bestScore = perfData.valid.perf;
                bestStrf = strf;
                bestIter = k;
            end

            clear perfData;

        end
        clear modelParamsCopy;
        fprintf('Found best STRF at iteration %d\n', bestIter);

        %% Set the best STRF as the trained strf
        modelParamsTrained.w1 = bestStrf;
    else       
        bestIter = -1;
        bestStrf = [];
        bestScore = -1;        
    end
    
    
    %% compute final performance with best STRF
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
        