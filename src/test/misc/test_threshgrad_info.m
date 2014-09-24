function test_threshgrad_info()

    %cellNames = {'yy1617_4_A', 'pupu2122_2_B'};
    cellNames = {'pupi0414_10', 'gg0116_1_A', 'obla1305_7', 'yy1617_5_B'};
    preprocNames = {'stft.nstd_6.fband_125', 'lyons','surprise.dfw_3.dtw_3.dg_4.Con'};
    
    for k = 1:length(cellNames)
        for m = 1:length(preprocNames)
            
            pname = preprocNames{m};
            if strcmp(pname, 'lyons')
                [earQ, step] = get_best_lyons_params(cellNames{k}, 'Con');
                pname = sprintf('lyons.agc_1.earQ_%d.step_%0.2f', earQ, step);
            end
            
            fprintf('Cell: %s  |  preproc: %s\n', cellNames{k}, pname);
            test_threshgrad_info_single(cellNames{k}, pname, 0.50);
            test_threshgrad_info_single(cellNames{k}, pname, 0.85);          
        end
    end

end


function test_threshgrad_info_single(cellName, preprocName, threshold)

    global globDat;

    h5 = h5utils();
    
    %cellName = 'yy1617_4_A';
    cellDir= fullfile('/auto/k6/mschachter/pystrfs/units', cellName);
    unitFile = fullfile(cellDir, 'unit.h5');
    outputDir = fullfile(cellDir, 'output');    
    
    %preprocName = 'stft.nstd_6.fband_125';
    preprocDir = '/auto/k6/mschachter/pystrfs/preproc';
    preprocFile = fullfile(preprocDir, sprintf('%s.h5', preprocName));
    
    stimClass = 'Con';
    
    cvRuns = 0:9;
    
    allGroups = 1:20;
    cvGroups = cv_partition(10, allGroups);
        
    for k = 1:length(cvRuns)
        
        
        %% get stimulus response data from file
        cvg = cvGroups{k};
        trainingGroups = cvg.training;
        validationGroups = cvg.validation;
        
        stimRespData = get_stimresp_data(unitFile, preprocFile, stimClass);
        
        %% get optimization data from file
        fname = sprintf('threshgrad.Con.%s.%0.2f.%d.h5', preprocName, threshold, cvRuns(k));        
        outputFile = fullfile(outputDir, fname);
        fprintf('\t%s\n', outputFile);
        fid = h5.open(outputFile);
    
        groupIndex = stimRespData.groupIndex;
        holdOuts = h5.get_ds(fid, '/data/holdout_subgroups');
        trainingIndex = h5.get_ds(fid, '/data/training_index');
        validationIndex = h5.get_ds(fid, '/data/validation_index');
        
        trStrs = h5.get_attr(fid, '/model', 'transforms');
        transforms = regexp(strtrim(trStrs), ',', 'split');
        
        strfs = h5.get_ds(fid, '/model/cv_weights');
        biases = h5.get_ds(fid, '/model/cv_bias');
        
        lastStrf = h5.get_ds(fid, '/model/weights');
        lastBias = h5.get_ds(fid, '/model/bias');
        
        numIters = h5.get_ds(fid, '/opt/num_iters');
        
        h5.close(fid);
        
        %% transform stimulus
        wholeStim = stimRespData.wholeStim;
        for m = 1:length(transforms)  
            tf = transforms{m};
            if ~isempty(tf)
                tfunc = sprintf('transform_%s(wholeStim)', tf);
                wholeStim = eval(tfunc);        
            end
        end
        stimRespData.wholeStim = wholeStim;
        
        
        %% Initialize strflab global variables with our stim and responses        
        strfData(stimRespData.wholeStim, stimRespData.wholeResp, stimRespData.groupIndex);

        %% Initialize a linear model
        strfLength = size(lastStrf, 2);
        strfDelays = 0:(strfLength-1);
        modelParams = linInit(stimRespData.numChannels, strfDelays, 'linear');
        
        meanStrf = squeeze(mean(strfs, 1));
        stdStrf = squeeze(std(strfs, 1));
        meanBias = mean(biases);
        
        infoFreqCutoff = 90;
        infoWindowSize = 0.200;    
        
            
        modelParams = linInit(stimRespData.numChannels, strfDelays, 'linear');            
        modelParams.w1 = lastStrf;
        modelParams.b1 = lastBias;

        [modelResp, rawModelResp] = compute_response(stimRespData, modelParams);
        
        fid = h5.open(outputFile);
        h5.set_ds(fid, '/model', 'response', modelResp);
        h5.set_ds(fid, '/model', 'response_raw', rawModelResp);
        h5.close(fid);        
        
        compute_response_info_coherence(outputFile);        
        
    end
    
end
   