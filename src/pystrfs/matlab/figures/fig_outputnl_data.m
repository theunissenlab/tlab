function fig_outputnl_data(unitName, stimClass, preprocType)

    switch preprocType
        case 'stft'
            preprocDesc = 'stft.nstd_6.fband_125';
        case 'lyons'
            preprocDesc = 'lyons.best';
        case 'surprise'
            preprocDesc = sprintf('surprise.dfw_3.dtw_3.dg_4.%s', stimClass);
    end
    
    unitDir = fullfile('/auto/k6/mschachter/pystrfs/units', unitName);
    outputDir = fullfile(unitDir, 'output');
    
    responseFileBase = sprintf('directfit.%s.%s.h5', stimClass, preprocDesc);
    responseFileDists = sprintf('directfit.nl_dists.%s.%s.h5', stimClass, preprocDesc);
    responseFileSpline = sprintf('directfit.nl_spline.%s.%s.h5', stimClass, preprocDesc);
    
    mdata = get_model_data(fullfile(outputDir, responseFileBase));
    mdata_dists = get_model_data(fullfile(outputDir, responseFileDists));
    mdata_spline = get_model_data(fullfile(outputDir, responseFileSpline));
    
    mdatas_desc = {'base', 'dists', 'spline'};
    mdatas = {mdata, mdata_dists, mdata_spline};
    
    minX = min(mdata.model.rawResponse);
    maxX = max(mdata.model.rawResponse);
    x = linspace(minX, maxX, 150);
    
    for k = 1:length(mdatas)
        m = mdatas{k};
        infoBound = m.coherence.bound.info_mean;
        infoTrain = m.coherence.training.info_mean;
        infoValid = m.coherence.validation.info_mean;
        
        perfTrain = infoTrain / infoBound;
        perfValid = infoValid / infoBound;
        
        fprintf('Model Info: %s\n', mdatas_desc{k});
        fprintf('\tUpper Bound: %0.1f bits\n', infoBound);
        fprintf('\tTraining: %0.1f bits, ratio=%0.3f\n', infoTrain, perfTrain);
        fprintf('\tValidation: %0.1f bits, ratio=%0.3f\n', infoValid, perfValid);
        fprintf('\n');
    end
    
    s = mdata.model.strf;
    cs = max(abs(s(:)));
    
    figure; hold on;
    subplot(2, 1, 1); hold on;
    imagesc(mdata.model.strf); axis tight;
    caxis([-cs cs]);
    title(sprintf('%s | %s | %s', unitName, stimClass, preprocType));
    
    subplot(2, 1, 2); hold on;
    plot(x, fnval(mdata_dists.model.outputNL, x), 'b-');
    plot(x, fnval(mdata_spline.model.outputNL, x), 'k-');
    legend('Dist', 'Spline');
    axis([minX maxX 0 2]);
    
        
    distinfo = mdata_dists.model.nlinfo;
    x = distinfo.x;    
    px_mean = mean(distinfo.px, 1);
    px_std = std(distinfo.px, 1);    
    pxspike_mean = mean(distinfo.pxspike, 1);
    pxspike_std = std(distinfo.pxspike, 1);    
    pxnospike_mean = mean(distinfo.pxnospike, 1);
    pxnospike_std = std(distinfo.pxnospike, 1);    
    pspikex_mean = mean(distinfo.pspikex, 1);
    pspikex_std = std(distinfo.pspikex, 1);    
    pnospikex_mean = mean(distinfo.pnospikex, 1);
    pnospikex_std = std(distinfo.pnospikex, 1);
    
    figure; hold on;    
    subplot(3, 1, 1); hold on;
    hist(mdata.model.rawResponse, 150);
    
    subplot(3, 1, 2); hold on;
    errorbar(x, px_mean, px_std, 'k-');
    errorbar(x, pxspike_mean, pxspike_std, 'r-');
    errorbar(x, pxnospike_mean, pxnospike_std, 'b-');
    legend('p(x)', 'p(x|spike)', 'p(x|nospike)');
    axis tight;
    
    subplot(3, 1, 3); hold on;
    errorbar(x, pspikex_mean, pspikex_std, 'r-');
    errorbar(x, pnospikex_mean, pnospikex_std, 'b-');
    axis([min(x) max(x) 0 2]);
    legend('p(spike|x)', 'p(nospike|x)');
    