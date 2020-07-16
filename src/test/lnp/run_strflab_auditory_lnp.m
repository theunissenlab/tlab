function outputFilePath = run_strflab_auditory_lnp(dataDir, strfLength, preprocType, nlType, optType, optOptions, desc, srFiles)

    global globDat

    %% set default parameters
    if nargin < 2
        strfLength = 40;
    end
    if nargin < 3
        preprocType = 'ft';
    end
    if nargin < 4
        nlType = 'exponential';
    end
    if nargin < 5
        optType = 'GradDesc';
    end
    if nargin < 6
        optOptions = struct;
        optOptions.display = 1;
        optOptions.coorDesc = 0;
        optOptions.earlyStop = 0;
        optOptions.stepSize = 1e-4;
        optOptions.nDispTruncate = 0;
        optOptions.maxIter = 10000;
        optOptions.lineSearch = 0;
    end
    if nargin < 7
        desc = 'default';
    end
    if nargin < 8
        srFiles = -1;
    end
    
    tic;
    
    
    %% check to see if the optimization has already been run
    outputDir = fullfile(dataDir, 'output');
    [s,mess,messid] = mkdir(outputDir);
    outputFileName = sprintf('strflab.lnp.%s.%s.%s.%s.mat', preprocType, nlType, optType, desc);
    outputFilePath = fullfile(outputDir, outputFileName); 
    fid = fopen(outputFilePath);
    if fid ~= -1
        fclose(fid);
        fprintf('Output file already exists for parameters, returning...\n');
        return;
    end
    
    %% get list of stim and response files
    if isstruct(srFiles)
        stimFiles = srFiles.stimFiles;
        respFiles = srFiles.respFiles;
    else
        stimFiles = get_filenames(dataDir, 'stim[0-9]', 1);    
        respFiles = get_filenames(dataDir, 'spike[0-9]*', 1);
    end
    
    %% preprocess stimuli and responses
    preprocDir = fullfile(dataDir, 'preproc');
    [s,mess,messid] = mkdir(preprocDir);  

    stimOutPattern = ['stim.preproc-' preprocType '.%d.mat'];
    respOutPattern = ['resp.preproc-' preprocType '.%d.mat'];
    
    srData = preprocess_sound(stimFiles, respFiles, preprocType, struct, preprocDir, stimOutPattern, respOutPattern);

    %% convert stim/response representation to strflab representation
    [allstim, allresp, groupIndex] = srdata2strflab_lnp(srData);
    
    %% subtract mean from stim
    stimNormalize = 0;
    if stimNormalize
        fprintf('Subtracting off mean stimuli and scaling by std dev...\n');
        numTimePoints = size(allstim, 1);
        stimStd = std(allstim);
        for k = 1:numTimePoints
            allstim(k, :) = allstim(k, :) - srData.stimAvg';
            allstim(k, :) = allstim(k, :) ./ stimStd;
        end        
    end
    
    usePsth = 1;
    if ~usePsth
        fprintf('Flattening PSTH...\n');
        allresp(allresp > 0) = 1;
    end

    %% Put stimulus and response and group assignments into globDat
    strfData(allstim, allresp, groupIndex);

    %% Set options for optimization method
    trnMethod = ['trn' optType];
    options = eval(trnMethod);
    options.display = 1;

    %override default options with user supplied options
    optFieldNames = fieldnames(optOptions);
    for k = 1:length(optFieldNames)
        fname = optFieldNames{k};
        options.(fname) = optOptions.(fname);
    end

    %% Initialize linear model
    modelParams = lnpInit(size(allstim,2), 0:strfLength);
    
    modelParams.regularize = 1;
    modelParams.regularizeWeight = 10;

    %% set the output nonlinearity
    modelParams.outputNL = nlType;

    % train the strf
    [modelParamsTrained, options] = strfOpt(modelParams, [], options);
    
    %% try prediction using the held out stimulus
    [modelParamsTrained, predResp] = strfFwd(modelParamsTrained, []);
    %scale response
    realResp = allresp;
    predResp = (predResp / max(predResp))*max(realResp);
    
    %compute coherence for validation
    sampleRate = 1000;
    cStruct = compute_coherence_mean(predResp, realResp', sampleRate);
    
    elapsedTime = toc;

    save(outputFilePath, 'cStruct', 'modelParamsTrained', 'options', 'preprocType', 'nlType', 'desc', 'elapsedTime');
