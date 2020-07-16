function test_threshgrad_filters()

    global globDat;

    h5 = h5utils();
    
    cellName = 'yy1617_4_A';
    cellDir= fullfile('/auto/k6/mschachter/pystrfs/units', cellName);
    unitFile = fullfile(cellDir, 'unit.h5');
    outputDir = fullfile(cellDir, 'output');    
    
    preprocName = 'stft.nstd_6.fband_125';
    preprocDir = '/auto/k6/mschachter/pystrfs/preproc';
    preprocFile = fullfile(preprocDir, sprintf('%s.h5', preprocName));
    
    stimClass = 'Con';
    
    %cvRuns = 0:9;
    cvRuns = [0];
    
    allGroups = 1:20;
    cvGroups = cv_partition(10, allGroups);
        
    for k = 1:length(cvRuns)
        
        
        %% get stimulus response data from file
        cvg = cvGroups{k};
        trainingGroups = cvg.training;
        validationGroups = cvg.validation;
        
        stimRespData = get_stimresp_data(unitFile, preprocFile, stimClass);
        
        %% get optimization data from file
        %fname = sprintf('threshgrad.Con.stft.nstd_6.fband_125.0.50.%d.h5', cvRuns(k));
        fname = sprintf('threshgrad.Con.stft.nstd_6.fband_125.0.85.%d.h5', cvRuns(k));
        outputFile = fullfile(outputDir, fname);
        fid = h5.open(outputFile);
    
        groupIndex = stimRespData.groupIndex;
        holdOuts = h5.get_ds(fid, '/data/holdout_subgroups')
        
        trStrs = h5.get_attr(fid, '/model', 'transforms');
        transforms = regexp(trStrs, ',', 'split');
        
        strfs = h5.get_ds(fid, '/model/cv_weights');
        biases = h5.get_ds(fid, '/model/cv_bias');
        
        lastStrf = h5.get_ds(fid, '/model/weights');
        lastBias = h5.get_ds(fid, '/model/bias');
        
        numIters = h5.get_ds(fid, '/opt/num_iters');
        
        
        %% transform stimulus
        wholeStim = stimRespData.wholeStim;
        for m = 1:length(transforms)       
            fprintf('Performing %s transform...\n', transforms{m});
            tfunc = sprintf('transform_%s(wholeStim)', transforms{m});
            wholeStim = eval(tfunc);        
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
        
        stds = [0 0.25 0.5 0.75 1 1.5 2];        
        cvPerfs = zeros(size(stds));        
        
        %% compute information for each std value
        for m = 1:length(stds)
            
            modelParamsFilt = linInit(stimRespData.numChannels, strfDelays, 'linear');            
            filteredStrf = df_fast_filter_filter(meanStrf, stdStrf, stds(m));
            modelParamsFilt.w1 = filteredStrf;
            modelParamsFilt.b1 = meanBias;
            
            [modelParamsFilt, predResp] = linFwd(modelParamsFilt, 1:length(groupIndex));
            predResp(isnan(predResp)) = 0;
            predResp(predResp < 0) = 0;
            
            psth = stimRespData.wholeResp;
            
            figure; hold on;
            subplot(2, 1, 1); hold on;
            imagesc(meanStrf); axis tight;
            cval = max(abs(meanStrf(:)));
            caxis([-cval cval]);
            title('Mean STRF');
            subplot(2, 1, 2); hold on;
            imagesc(filteredStrf); hold on; axis tight;
            cval = max(abs(filteredStrf(:)));
            caxis([-cval cval]);
            title(sprintf('Filtered STRF: std=%d', stds(m)));
            
            figure; hold on;
            subplot(2, 1, 1); hold on;
            plot(psth); axis tight;
            subplot(2, 1, 2); hold on;
            plot(predResp); axis tight;            
            
            cStruct = compute_coherence_mean(predResp, psth, stimRespData.sampleRate, infoFreqCutoff, infoWindowSize);
            cvPerfs(m) = cStruct.info;        
            fprintf('\tstd=%d, info=%0.2f bits\n', stds(m), cStruct.info);
            
            %{
            nPrts = size(holdOuts, 1);
            lens = zeros(nPrts, 1);
            preds = cell(nPrts, 1);
            psths = cell(nPrts, 1);           
            
            for p = 1:nPrts
               
                %% get training and holdout indicies
                vIndx = findIdx(holdOuts(p, :), groupIndex);
                gindx = groupIndex;
                gindx(vIndx) = 0;
                tIndx = find(gindx > 0);

                modelParamsFilt = linInit(stimRespData.numChannels, strfDelays, 'linear');
                strfP = squeeze(strfs(p, :, :));
                filteredStrf = df_fast_filter_filter(strfP, stdStrf, stds(m));
                
                modelParamsFilt.w1 = filteredStrf;
                modelParamsFilt.b1 = biases(p);
                
                lens(p) = length(vIndx);
                psths{p} = stimRespData.wholeResp(vIndx);
                [modelParamsFilt, predResp] = linFwd(modelParamsFilt, vIndx);
                predResp(isnan(predResp)) = 0;
                predResp(predResp < 0) = 0;
                preds{p} = predResp;
                
            end
            
            %% concatenate responses
            concatPreds = zeros(sum(lens), 1);
            concatPsth = zeros(sum(lens), 1);
            lastIndex = 0;
            for p = 1:nPrts
                si = lastIndex + 1;
                ei = lastIndex + lens(p);
                concatPreds(si:ei) = preds{p};
                concatPsth(si:ei) = psths{p};            
                lastIndex = lastIndex + lens(p);
            end
            
            
            %% compute info
            cStruct = compute_coherence_mean(concatPreds, concatPsth, stimRespData.sampleRate, infoFreqCutoff, infoWindowSize);
            cvPerfs(m) = cStruct.info;        
            fprintf('\tstd=%d, info=%0.2f bits\n', stds(m), cStruct.info);
            %}
                        
        end
        
        
    end
    
end
   