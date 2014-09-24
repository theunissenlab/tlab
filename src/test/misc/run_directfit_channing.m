function run_directfit_channing()

    dataDir = '/auto/fdata/mschachter/data/channing/clean';
    
    addpath('/auto/fhome/channing/workspace/databrowser/browser/model/analysis/strflab');    
    [audioFiles, allSpikeTrials, represp_ids] = strflab_get_response_data(307564);
        
    %% get and compare stimuli from strfpak and strflab
    wavFiles = cell(length(audioFiles), 1);
    spakStimFiles = cell(length(audioFiles), 1);
    spakSpikeFiles = cell(length(audioFiles), 1);
    
    for k = 1:length(audioFiles)
        [pathstr, name, ext] = fileparts(audioFiles{k});
        if ~strcmp('bos', name)
            wavFile =  sprintf('%s.wav', name);
            spakFile = sprintf('%s.mat', name);            
            sindx = str2num(name(5:end));            
        else
            wavFile = 'stim45.wav';
            spakFile = 'stim45.mat';
            sindx = 45;
        end
        spakStimFiles{k} = fullfile(dataDir, spakFile);
        wavFiles{k} = fullfile(dataDir, wavFile);
        spakSpikeFile = sprintf('Spike%d.mat', sindx);
        spakSpikeFiles{k} = fullfile(dataDir, spakSpikeFile);
        %fprintf('%d %s | %s | %s | %s\n', k, spakFile, wavFiles{k}, audioFiles{k}, spakSpikeFile);
    end

    preprocStimParams = struct;      %create preprocessing param structure
    preprocStimParams.tfType = 'stft'; %use short-time FT
    tfParams = struct;               %create time-frequency params
    tfParams.high_freq = 8000;       %specify max freq to analyze
    tfParams.low_freq = 250;         %specify min freq to analyze
    tfParams.log = 1;                %take log of spectrogram
    tfParams.dbnoise = 80;           %cutoff in dB for log spectrogram, ignore anything below this
    tfParams.refpow = 0;             %reference power for log spectrogram, set to zero for max of spectrograms across stimuli
    tfParams.overwrite = 1;
    preprocStimParams.tfParams = tfParams;

    preprocStimParams.outputDir = fullfile(dataDir, 'preproc');

    [wholeStim, groupIndex, stimInfo, preprocStimParams] = preprocSound(wavFiles, preprocStimParams);
    fprintf('z-scoring stimulus...\n');    
    wholeStim = norm_std_mean(wholeStim);
    

    %{
    wholeStim = [];
    stimLengths = zeros(length(spakStimFiles), 1);
    groupIndex = [];
    for k = 1:length(spakStimFiles)
       
        sfile = spakStimFiles{k};
        vars = load(sfile);
        stim = vars.outSpectrum;
        
        slen = size(stim, 2);
        stimLengths(k) = slen/1000;
        wholeStim = [wholeStim; stim'];
        gindx = ones(1, slen)*k;
        groupIndex = [groupIndex gindx];
    end
    
    stimInfo.stimLengths = stimLengths;
    stimInfo.numStimFeatures = size(wholeStim, 2);
    
    size(wholeStim)
    size(groupIndex)
    %}
    
    groups = unique(groupIndex);
    %{
    for k = groups
       
        spakVars = load(spakStimFiles{k});        
        spakStim = spakVars.outSpectrum;
        clear spakVars;
        
        slabStim = wholeStim(groupIndex == k, :)';
        
        slabLen = size(slabStim, 2);
        spakLen = size(spakStim, 2);
        
        lendiff = slabLen - spakLen;
        
        figure; hold on;
        subplot(2, 1, 1); hold on;
        imagesc(slabStim); axis tight;
        colorbar;
        title(wavFiles{k});
        
        subplot(2, 1, 2); hold on;
        imagesc(spakStim); axis tight;
        colorbar;
        title(spakStimFiles{k});
                
        fprintf('slab=%d | spak=%d | diff=%d\n', slabLen, spakLen, lendiff);
    end
    %}
    
    %% get responses from strfpak    
    allSpikeTrials2 = cell(length(spakSpikeFiles), 1);
    for k = 1:length(spakSpikeFiles)    
        spikeVars = load(spakSpikeFiles{k});
        allSpikeTrials2{k} = spikeVars.rawResp;
    end
    
    
    %% create PSTH
    preprocRespParams = struct;         
    preprocRespParams.units = 'ms';     
    preprocRespParams.binSize = 0.001;  
    
    preprocRespParams.stimLengths = stimInfo.stimLengths;
    [wholeResponse, groupIndex, respInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);
    %[wholeResponse2, groupIndex2, respInfo2, preprocRespParams] = preprocSpikeResponses(allSpikeTrials2, preprocRespParams);
    
%{    
    for k = groups
       
        resp1 = wholeResponse(groupIndex == k);
        resp2 = wholeResponse2(groupIndex2 == k);
        
        len = min(length(resp1), length(resp2));
        len1 = length(resp1);
        len2 = length(resp2);
        slen = round(stimInfo.stimLengths(k)*1000);
        fprintf('stimLen=%d | len1=%d | len2=%d\n', slen, len1, len2);
        
        figure; hold on;
        plot(resp1(1:len), 'k-', 'LineWidth', 2); axis tight;
        plot(resp2(1:len), 'r-'); axis tight;
        title(spakSpikeFiles{k});
        
    end
  %}  
    
    %% get training and early stopping dataset indicies
    trainingIndex = findIdx(1:45, groupIndex);
    %validationIndex = findIdx(validationGroups, groupIndex);
    
    
    %% Initialize strflab global variables with our stim and responses
    global globDat
    strfData(wholeStim, wholeResponse, groupIndex);


    %% Initialize a linear model
    strfLength = 100;
    strfDelays = 0:(strfLength-1);
    modelParams = linInit(stimInfo.numStimFeatures, strfDelays, 'linear');
    modelParams.b1 = mean(wholeResponse);
    
    
    %% Initialize Direct Fit options
    optOptions = trnDirectFit();
    %optOptions.tolerances = [0.500 0.100 0.050 0.010 0.005 1e-03 5e-04 1e-04 5e-05];
    optOptions.tolerances = [0.1 0.05 0.010 0.005 1e-03 5e-04 1e-04 5e-05];
    %optOptions.sparsenesses = [0 1 2 4 6 8 10 12];
    optOptions.sparsenesses = [0 1 2 6];
    optOptions.infoFreqCutoff = 90;
    optOptions.infoWindowSize = 0.250;
       
    %% run direct fit
    [modelParamsTrained, optOptions] = strfOpt(modelParams, trainingIndex, optOptions);
    
    
    save(fullfile(dataDir, 'strflab.output.mat'), 'modelParamsTrained');
    
    figure; hold on;
    imagesc(modelParamsTrained.w1); axis tight;
