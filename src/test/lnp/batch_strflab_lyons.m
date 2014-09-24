function batch_strflab_lyons()

    rootdir = '/auto/fdata/mschachter/data';
    stimdir = '/auto/fdata/mschachter/data/all_stims';

    if nargin < 3
        %sdNames = {'yy1617_4_A'}; %shitty
        %sdNames = {'gg0116/5_A'}; %shitty
        sdNames = {'pupu1203_4_A'};
        %sdNames = {'oo0909_2_A'};
        %sdNames = {'pipi1112_1_A'};        
    end
    
    dirNames = cell(length(sdNames), 1);
    for k = 1:length(sdNames)
        dirNames{k} = fullfile(rootdir, sdNames{k});
    end
        
    strfLength = 59;
        
    datasets = find_datasets(dirNames{k}, stimdir, 'conspecific|zfsongs');
    srFiles = datasets{1}.srPairs;
    dataDir = datasets{1}.dirname;

    holdOutIndex = -1;

    preprocType = 'lyons';
    
    %stepParams = [0.75 0.5 0.25];
    stepParams = [0.75];
    earQParams = [8];
    %earQParams = [14 16 18];
    %stepParams = [0.5];
    %earQParams = [8];
    
    paramStructs = cell(length(stepParams)*length(earQParams), 1);
    indx = 1;
    for k = 1:length(stepParams)
        for j = 1:length(earQParams)

            step = stepParams(k);
            earQ = earQParams(j);
            pname = ['lyons.step' num2str(step) '.earQ' num2str(earQ)];
            
            lyonsParams = struct;
            lyonsParams.name = pname;
            lyonsParams.low_freq = 250;
            lyonsParams.high_freq = 9000;
            lyonsParams.earQ = earQ;
            lyonsParams.agc = 1;
            lyonsParams.differ = 1;
            lyonsParams.tau = 3;
            lyonsParams.step = step;
            
            paramStructs{indx} = lyonsParams;
            indx = indx + 1;
            
        end
    end
    
    
    preprocType = 'ft';
    ftparam = struct;
    ftparam.name = 'ft';
    paramStructs = {ftparam};
    
        
    %% preprocess stimuli and responses
    stimFiles = srFiles.stimFiles;
    respFiles = srFiles.respFiles;

    preprocDir = fullfile(dataDir, 'preproc');
    [s,mess,messid] = mkdir(preprocDir);  
    
    for k = 1:length(paramStructs)

        preprocParams = paramStructs{k};        
        pname = preprocParams.name;
        stimOutPattern = ['stim.preproc-' pname '.%d.mat'];
        respOutPattern = ['resp.preproc-' pname '.%d.mat'];

        srData = preprocess_sound(stimFiles, respFiles, preprocType, preprocParams, preprocDir, stimOutPattern, respOutPattern);
        nStimChannels = srData.nStimChannels;

        %% Linear Model + direct fit            
        strfDelays = -15:strfLength;
        linModelParams = linInit(nStimChannels, strfDelays);
        
        preprocOptions = struct;
        
        optOptions = struct;
        optOptions.optName = 'DirectFit';

        desc = pname;
        overwrite = 1;
        outputFilePath = run_strflab_auditory(dataDir, srData, preprocOptions, linModelParams, optOptions, desc, holdOutIndex, overwrite);
        strflab_auditory_compute_predictions(outputFilePath, overwrite);
    end
    
    