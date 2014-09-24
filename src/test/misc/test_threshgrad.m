function test_threshgrad()
    
    h5 = h5utils();
    
    cellName = 'yy1617_4_A';
    cellDir= fullfile('/auto/k6/mschachter/pystrfs/units', cellName);
    outputDir = fullfile(cellDir, 'output');    
    
    %cvRuns = 0:9;
    cvRuns = [0];
    
    for k = 1:length(cvRuns)
        fname = sprintf('threshgrad.Con.stft.nstd_6.fband_125.0.50.%d.h5', cvRuns(k));
        outputFile = fullfile(outputDir, fname);
        fid = h5.open(outputFile);
    
        strfs = h5.get_ds(fid, '/model/cv_weights');
        biases = h5.get_ds(fid, '/model/cv_bias');
        
        lastStrf = h5.get_ds(fid, '/model/weights');
        lastBias = h5.get_ds(fid, '/model/bias');
        
        trainingIndex = h5.get_ds(fid, '/data/training_index');
        validationIndex = h5.get_ds(fid, '/data/validation_index');
        
        numIters = h5.get_ds(fid, '/opt/num_iters');
        
        numStrfs = size(strfs, 1);
        figure('Name', sprintf('CV Run %d', k)); hold on;
        for m = 1:numStrfs           
            strf = squeeze(strfs(m, :, :));            
            bias = biases(m, :, :);
            c = max(abs(strf(:)));
            
            subplot(numStrfs, 1, m); hold on;
            imagesc(strf); axis tight;
            caxis([-c c]);            
            title(sprintf('# iter=%d, Bias=%f', numIters(m), bias));
        end
        
        strfMean = squeeze(mean(strfs, 1));
        strfStd = squeeze(std(strfs, 1));
        
        
        
        
        
        
        figure('Name', sprintf('CV Run %d', k)); hold on;
        
        subplot(3, 1, 1); hold on;
        imagesc(strfMean); axis tight;
        c = max(abs(strfMean(:)));
        caxis([-c c]);
        title(sprintf('STRF Mean: bias=%f +/- %f', mean(biases), std(biases)));
        
        subplot(3, 1, 2); hold on;
        imagesc(strfStd); axis tight;
        c = max(abs(strfStd(:)));
        caxis([-c c]);
        title('STRF std');
        
        subplot(3, 1, 3); hold on;
        imagesc(lastStrf); axis tight;
        c = max(abs(lastStrf(:)));
        caxis([-c c]);
        title(sprintf('Last STRF, bias=%f', lastBias));
        
        [strfMeanFiltered, strfMeanFilter] = df_fast_filter_filter(strfMean, strfStd, 1);
        [strfFiltered, strfFilter] = df_fast_filter_filter(lastStrf, strfStd, 1);
        figure('Name', sprintf('CV Run %d', k)); hold on;
        
        subplot(4, 1, 1); hold on;
        imagesc(strfMeanFilter); axis tight;
        c = max(abs(strfMeanFilter(:)));
        caxis([-c c]);
        colorbar;
        title('Mean Filter');
        
        subplot(4, 1, 2); hold on;
        imagesc(strfMeanFiltered); axis tight;
        c = max(abs(strfMeanFiltered(:)));
        caxis([-c c]);
        colorbar;
        title('Mean Filtered STRF');
        
        subplot(4, 1, 3); hold on;
        imagesc(strfFilter); axis tight;
        c = max(abs(strfFilter(:)));
        caxis([-c c]);
        colorbar;
        title('Filter');
        
        subplot(4, 1, 4); hold on;
        imagesc(strfFiltered); axis tight;
        c = max(abs(strfFiltered(:)));
        caxis([-c c]);
        colorbar;
        title('Filtered STRF');        
        
    end
    