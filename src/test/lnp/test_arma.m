function test_arma()

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

    %% fit AR model
    filterLength = 50;
    m = ar(allresp, filterLength);
    
    srFilter = get(m, 'a');
    srFilter = srFilter(2:end);
    
    %figure; hold on;
    %plot(srFilter); axis tight;
    
    predResp = zeros(size(allresp));
    for k = 1:length(predResp)
      
      
      
    end
    
    predResp = conv(srFilter, allresp);
    figure; hold on;
    plot(allresp, 'k-');
    plot(predResp, 'b-');
    axis tight;
    
    
    