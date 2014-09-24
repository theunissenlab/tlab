function test_lnpsr2_real()

    cellDir = '/auto/fdata/mschachter/data/pipi1112_1_A';
    stimsDir = '/auto/fdata/mschachter/data/all_stims';
    
    datasets = find_datasets(cellDir, stimsDir, 'conspecific|zfsongs');
    srFiles = datasets{1}.srPairs;
    dataDir = datasets{1}.dirname;

    preprocType = 'ft';

    %% preprocess stimuli and responses
    stimFiles = srFiles.stimFiles;
    respFiles = srFiles.respFiles;

    preprocDir = fullfile(dataDir, 'preproc');
    [s,mess,messid] = mkdir(preprocDir);  

    stimOutPattern = ['stim.preproc-' preprocType '.%d.mat'];
    respOutPattern = ['resp.preproc-' preprocType '.%d.mat'];

    srData = preprocess_sound(stimFiles, respFiles, preprocType, ...
			      struct, preprocDir, stimOutPattern, respOutPattern);
    numStimChannels = srData.nStimChannels;
    
    [allstim, allresp, groupIndex] = srdata2strflab(srData, 0);
    clear allresp; %not used
    
    %% load preexisting STRF
    strfFile = fullfile(dataDir, 'output/strf.pipi1112_1_A.mat');
    vars = load(strfFile);
    strf = vars.strf;
    %strf = strfGain*simple_gabor_strf(numStimChannels, stimDelays, 5, 30); 
    clear vars;
    numDelays = size(strf, 2);
    stimDelays = 0:(numDelays - 1);
    strfGain = 2e2;
    strf = strfGain*strf;
    
    %% create fake SR filter
    spikeDelays = [];
    
    %% construct LNP model to generate response
    ignoreBias = 1;
    nlType = 'exponential';    
    sampleRate = 1000;
    modelParamsReal = lnpsr2Init(numStimChannels, stimDelays, nlType, spikeDelays, sampleRate, ignoreBias);
    
    modelParamsReal.w1 = strf;
    
    %% create response from fixed LNP model w/ real stimuli
    global globDat;
    nTimePoints = size(allstim, 1);
    strfData(allstim, zeros(1, nTimePoints), groupIndex);
    
    numTrials = 20;
    datIdx = 1:nTimePoints;
    [modelParamsReal, modelResponseReal, fullResponse] = lnpsr2Fwd(modelParamsReal, datIdx, numTrials);
        
    plotRealResps = 0;
    if plotRealResps      
      figure; hold on;
      imagesc(strf); axis tight; colorbar;
      title('STRF');
      
      stimCurMax = max(fullResponse.stimCurrent);
      stimCurMin = min(fullResponse.stimCurrent);
      for k = 1:length(srData.datasets)	
        ds = srData.datasets{k};
        tfrep = ds.stim.tfrep;        
        indx = find(groupIndex == k);
        stim = allstim(indx, :);
        stimCur = fullResponse.stimCurrent(indx);
        nonlinResp = fullResponse.nonlinearResponse(1, indx);
        nlAvg = mean(nonlinResp(:));
        nlStd = std(nonlinResp(:));
        nlMin = 1;
        nlMax = nlAvg + 4*nlStd;
        resp = modelResponseReal(indx);
        
        figure; hold on;
        
        subplot(4, 1, 1);
        %imagesc(stim'); axis tight;
        plot_tfrep(tfrep);
        title(sprintf('Stim %d', k));
        
        subplot(4, 1, 2);
        plot(stimCur);
        axis([1, length(stimCur), stimCurMin, stimCurMax]);
        
        subplot(4, 1, 3); hold on;
        plot(nonlinResp');
        axis([1, length(nonlinResp), nlMin, nlMax]);
        
        subplot(4, 1, 4);
        plot(resp, 'k-');	
        axis([1, length(resp), 0, 1]);
        title(sprintf('avg rate=%f', mean(resp)));
        
      end      
    end
    
    plotTrials = 0;
    if plotTrials
      
      indx = (groupIndex == k);
      fresp = fullResponse;
      for k = 1:fresp.numTrials
		
        figure; hold on;
        subplot(2, 1, 1); hold on;
        plot(fresp.stimCurrent(indx), 'r-');
        plot(fresp.spikeCurrent(indx), 'b-');
        axis tight;

        subplot(2, 1, 2); hold on;
        plot(fresp.nonlinearResponse(indx), 'k-');
        stimes = find(fresp.spikeTrials(k, indx) > 0) - 1;
        svals = ones(size(stimes))*0.8*max(fresp.nonlinearResponse(indx));
        plot(stimes, svals, 'ro');		
      end
      
      
    end

    spikeTimes = cell(numTrials, 1);
    for k = 1:numTrials
        st = fullResponse.spikeTrials(k, :);
        spikeTimes{k} = (find(st > 0) - 1) / sampleRate;
    end
    
    %% store spike times in global response for model construction
    strfData(allstim, spikeTimes, groupIndex);
    
    %% Initialize new LNP spike-response model to fit original
    nlType = 'exponential';
    modelParams = lnpsr2Init(numStimChannels, stimDelays, nlType, spikeDelays, sampleRate);
    
    modelParams.fixedStrf = 0;
    if modelParams.fixedStrf
        fprintf('Using fixed STRF...\n');
        modelParams.w1 = modelParamsReal.w1;
        modelParams.b1 = modelParamsReal.b1;
    end 
        
    usePerturbedInitialGuess = 0;
    if usePerturbedInitialGuess
        fprintf('Using perturbed initial guess...\n');
        rgain = 1e-3;
        strf = modelParamsReal.w1;
        modelParams.w1 = strf + randn(size(strf))*mean(mean(strf))*rgain;
        modelParams.b1 = modelParamsReal.b1 + randn()*modelParamsReal.b1*rgain;
        swts = modelParams.spikeResponseWeights;
        modelParams.spikeResponseWeights = swts + randn(size(swts))*mean(swts)*rgain;
    end
    
    useRandomInitialGuess = 1;
    if useRandomInitialGuess
        fprintf('Using random initial guess...\n');
        rgain = 1e-3;
        strf = modelParamsReal.w1;
        modelParams.w1 = randn(size(strf))*mean(mean(strf))*rgain;
        modelParams.b1 = randn()*modelParamsReal.b1*rgain;
        swts = modelParamsReal.spikeResponseWeights;
        modelParams.spikeResponseWeights = randn(size(swts))*mean(swts)*rgain;
    end
    
    modelParams.checkGrad = 0;    
    modelParams.useFDGrad = 0;
    
    %% run strflab optimization
    optOptions = trnSCG;
    optOptions.display = 1;
    optOptions.maxIter = 1000;
    [modelParamsTrained, optOptions] = strfOpt(modelParams, datIdx, optOptions);
       
    %% look at response
    [modelParamsTrained, predResp] = lnpsr2Fwd(modelParamsTrained, datIdx, numTrials);
    predResp = (predResp/max(predResp)) * max(modelResponseReal);
    
    realStrf = modelParamsReal.w1;
    predStrf = modelParamsTrained.w1;
    strfDiff = realStrf - predStrf;
    nstrfDiff = norm(strfDiff);
  
    realw0 = modelParamsReal.b1;
    predw0 = modelParamsTrained.b1;
    
    realSwts = modelParamsReal.spikeResponseWeights;
    predSwts = modelParamsTrained.spikeResponseWeights;
    
    figure; hold on;
    plot(modelResponseReal, 'k-', 'LineWidth', 2);
    plot(predResp, 'r-', 'LineWidth', 2);
    axis tight;
    legend('Real', 'Pred');
    
    figure; hold on;
    subplot(3, 1, 1);
    imagesc(realStrf); axis tight;
    title(sprintf('Real STRF, bias=%f', realw0));
    colorbar;
    
    subplot(3, 1, 2);
    imagesc(predStrf); axis tight;
    title(sprintf('Pred STRF, bias=%f', predw0));
    colorbar;
    
    subplot(3, 1, 3);
    imagesc(strfDiff); axis tight;
    title(sprintf('STRF diff, norm=%f', nstrfDiff));
    colorbar;  
    
    save(fullfile(cellDir, 'output'), ...
	 'modelParamsReal', 'modelParamsTrained', 'modelResponseReal', ...
	 'predResp');
%}    