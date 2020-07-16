function run_outputnl_sim(dataDir, cellName, preprocType, stimType, nlType, outputFilePath)

    %% set up directories and load STRFs
    cellDir = fullfile(dataDir, cellName);
    outputNLDir = fullfile(dataDir, 'outputNL');    
    strfsDir = fullfile(outputNLDir, 'strfs');
    
    strfFile = fullfile(strfsDir, sprintf('%s.mat', cellName));
    
    vars = load(strfFile);
    strfs = vars.strfs;    
    clear vars;
    
    %% run simulations for the preprocType,stimType pair and save outpu to file    
    strf = strfs.(preprocType).(stimType).strf;
    bias = strfs.(preprocType).(stimType).bias;            

    
    %% get stimulus from actual cell directory
    slabOutputDir = fullfile(cellDir, stimType, 'output');
    slabFileName = sprintf('strflab.tfType_%s.%s.mat', preprocType, stimType);
    slabFile = fullfile(slabOutputDir, slabFileName);
    vars = load(slabFile);            
    groupIndex = vars.preprocData.groupIndex;
    stim = vars.preprocData.stim;
    if strcmp('surprise', preprocType)
        specStim = vars.preprocData.stim;
        fprintf('Computing surprise of spectrogram...\n');
        surpriseParams = struct;
        surpriseParams.outputPath = fullfile(cellDir, stimType, 'preproc');
        surpriseParams.outputDesc = 'default';
        [surpriseStimLouder, surpriseStimQuieter, groupIndex, surpriseStimInfo, surpriseParams] = preprocSoundSurprise(specStim, groupIndex, surpriseParams);
        stim = [surpriseStimLouder; surpriseStimQuieter]';
    end    
    groups = unique(groupIndex);
    
    %% get mean spike rate
    meanResp = vars.preprocData.respInfo.mean;
    clear vars;
    
    
    %% pick training and validation data
    glen = length(groups);
    validationGroups = groups((glen-1):glen);
    trainingGroups = 1:(glen-2);
    

    %% set up GLPP model with specified output nonlinearity    
    sampleRate = 1000;
    numStimChannels = size(strf, 1);
    delays = 0:(size(strf, 2)-1);
    modelParams = glppInit(numStimChannels, delays, 'exponential', [], sampleRate, 0);
    
    modelParams.w1 = strf;
    modelParams.b1 = bias;

    nSamplePoints = size(stim, 1);
    resp = zeros(1, nSamplePoints);
    strfData(stim, resp, groupIndex);

    
    %% find ideal multiplier  to normalize strf so that response of model is equal to mean response
    strfMultiplier = find_multiplier(modelParams, meanResp);
    modelParams.w1 = strf*strfMultiplier;
    modelParams.b1 = bias*strfMultiplier;    
    
    %% simulate spike trials
    fprintf('Simulating spike trials...');
    numTrials = 20;
    
    [modelParams, modelResponse, fullResponse] = glppFwd(modelParams, 1:length(groupIndex), numTrials);

    spikeTrials = cell(numTrials, 1);
    for j = 1:numTrials       
        stimes = (find(fullResponse.spikeTrials(j, :))-1) / sampleRate;
        spikeTrials{j} = stimes;        
    end    
    allSpikeTrials = {spikeTrials};
    fprintf('done!\n');
    
    %% pick out training and validation sets
    trainingIndex = findIdx(trainingGroups, groupIndex);    
    validationIndex = findIdx(validationGroups, groupIndex);        
    
    linRespTrain = fullResponse.stimCurrent(trainingIndex);    
    psthTrain = modelResponse(trainingIndex);
    
    linRespValid = fullResponse.stimCurrent(validationIndex);
    psthValid = modelResponse(validationIndex);

    
    trainGroupIndex = groupIndex(trainingIndex);

    %% fit output nonlinearity on training data using 2 different methods
    tic;
    nlinfoDists = estimate_outputnl_from_dists(linRespTrain, psthTrain, numTrials, trainGroupIndex);    
    dtime = toc;
    fprintf('Estimation from spike-conditioned densities took %fs\n', dtime);
    
    tic;
    nlinfoSpline = estimate_outputnl_from_spline(linRespTrain, psthTrain, numTrials, trainGroupIndex);
    sptime = toc;
    fprintf('Estimation from spline took %fs\n', sptime);

    %% split PSTH for coherence measure
    infoFreqCutoff = 90;
    infoWindowSize = 0.250;
    respParams = struct;
    respParams.units = 's';
    respParams.split = 1;
    respParams.stimLengths = [nSamplePoints / sampleRate];
    [wholeResponse, groupIndex, respInfo, respParams] = preprocSpikeResponses(allSpikeTrials, respParams);

    
    %% compute reponses and performance
    fprintf('Computing responses...\n');
    nlnames = {'p(spike|x)', 'spline'};
    nlInfos = {nlinfoDists, nlinfoSpline};
    perfInfo = cell(length(nlInfos), 1);
    ofuncTypes = {'err', 'kl'};
    
    for k = 1:length(nlInfos)
        
        pinfo = struct;
        ni = nlInfos{k};
        
        pinfo.nlinfo = ni;
        pinfo.psthTrain = psthTrain;
        pinfo.psthValid = psthValid;
        pinfo.name = nlnames{k};        
        
        for m = 1:length(ofuncTypes)
            
            oftype = ofuncTypes{m};
            
            %% compute responses from p(spike|x) NL
            trainResp = fnval(ni.(oftype).outputNL, linRespTrain);
            validResp = fnval(ni.(oftype).outputNL, linRespValid);

            %% compute performances    
            [cBoundTrain, cTrain] = compute_coherence_full(trainResp, psthTrain, wholeResponse(1, trainingIndex), wholeResponse(2, trainingIndex), ...
                                                      sampleRate, numTrials, infoFreqCutoff, infoWindowSize);
            perfRatioTrain = cTrain.info / cBoundTrain.info;

            [cBoundValid, cValid] = compute_coherence_full(validResp, psthValid, wholeResponse(1, validationIndex), wholeResponse(2, validationIndex), ...
                                                      sampleRate, numTrials, infoFreqCutoff, infoWindowSize);
            perfRatioValid = cValid.info / cBoundValid.info;
            
            pinfo.(oftype).trainResp = trainResp;
            pinfo.(oftype).validResp = validResp;
            pinfo.(oftype).cBoundTrain = cBoundTrain;
            pinfo.(oftype).cBoundValid = cBoundValid;
            pinfo.(oftype).perfRatioTrain = perfRatioTrain;
            pinfo.(oftype).perfRatioValid = perfRatioValid;

        end
        
        perfInfo{k} = pinfo;        

    end
        
    %% save stuff to an output file
    outputFileName = sprintf('%s.%s.%s.%s.mat', cellName, preprocType, stimType, nlType);
    outputFile = fullfile(outputFilePath, outputFileName);

    save(outputFile, 'modelParams', 'numTrials', 'modelResponse', 'fullResponse', 'trainingIndex', 'validationIndex', ...
                     'linRespTrain', 'psthTrain', 'linRespValid', 'psthValid', 'perfInfo');

    fprintf('Done!\n');
    
    
    
end


%% bisection method to determine optimum STRF multiplier, i.e. one that
%% produces a response closest to targetMeanResp
function strfMultiplier = find_multiplier(modelParams, targetMeanResp)

    global globDat;

    numTrials = 10;
    origModelParams = modelParams;

    maxMult = 75;
    minMult = 0.5;
    
    a = minMult;
    b = maxMult;
    
    minInterval = 1;
    meanRespTargetDiff = 1e-3;
    
    converged = 0;
    while ~converged
        %% choose midpoint
        c = ((b - a) / 2) + a;
        
        %% use midpoint as multiplier
        mult = c;
        modelParams = origModelParams;
        modelParams.w1 = modelParams.w1*mult;
        modelParams.b1 = modelParams.b1*mult;
        
        [modelParams, modelResponse] = glppFwd(modelParams, 1:length(globDat.groupIdx), numTrials);    
        
        meanResp = mean(modelResponse);
        meanRespDiff = meanResp - targetMeanResp;
        fprintf('[%0.2f, %0.2f], c=%0.2f, meanResp=%0.4f, diff=%0.4f\n', a, b, c, meanResp, meanRespDiff);
        
        %% choose new interval
        if meanRespDiff >= 0
            b = c;
        else
            a = c;            
        end
            
        %% check for convergence
        converged = abs(meanRespDiff) < meanRespTargetDiff;
        converged = converged || (b-a) < minInterval;
    end
    
    strfMultiplier = mult;
    
end
    