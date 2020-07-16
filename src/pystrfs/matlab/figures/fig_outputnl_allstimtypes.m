function fig_outputnl_allstimtypes(allData)

    preprocTypes = {'stft', 'lyons', 'surprise'};   
    stimClasses = {'Con', 'Flatrip'};
    
    strfFig = figure('Name', sprintf('%s: STRFs', allData.unitName));
    distFig = figure('Name', sprintf('%s: Dists', allData.unitName));
    distFig2 = figure('Name', sprintf('%s: ST-Dists', allData.unitName));
    nlFig = figure('Name', sprintf('%s: Output NLs', allData.unitName));
    histFig = figure('Name', sprintf('%s: Hists', allData.unitName));
    
    linTrainPerfs = zeros(length(stimClasses), length(preprocTypes));
    linValidPerfs = zeros(length(stimClasses), length(preprocTypes));
    distTrainPerfs = zeros(length(stimClasses), length(preprocTypes));
    distValidPerfs = zeros(length(stimClasses), length(preprocTypes));
    splineTrainPerfs = zeros(length(stimClasses), length(preprocTypes));
    splineValidPerfs = zeros(length(stimClasses), length(preprocTypes));
    
    for k = 1:length(stimClasses)        
        stimClass = stimClasses{k};        
        
        for m = 1:length(preprocTypes)            
            
            spNum = (2*m - 1) + (k - 1);
            
            preprocType = preprocTypes{m};
            
            lin = allData.(stimClass).(preprocType).linear;
            dist = allData.(stimClass).(preprocType).dists;
            spline = allData.(stimClass).(preprocType).spline;
          
            linTrainPerfs(k, m) = lin.coherence.training.info_mean / lin.coherence.bound.info_mean;
            linValidPerfs(k, m) = lin.coherence.validation.info_mean / lin.coherence.bound.info_mean;
            distTrainPerfs(k, m) = dist.coherence.training.info_mean / dist.coherence.bound.info_mean;
            distValidPerfs(k, m) = dist.coherence.validation.info_mean / dist.coherence.bound.info_mean;
            splineTrainPerfs(k, m) = spline.coherence.training.info_mean / spline.coherence.bound.info_mean;
            splineValidPerfs(k, m) = spline.coherence.validation.info_mean / spline.coherence.bound.info_mean;
                        
            strf = lin.model.strf;
            ca = max(abs(strf(:)));
            
            distinfo = dist.model.nlinfo;
            x = distinfo.x;    
            minX = min(x); maxX = max(x);
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
            
            figure(strfFig);            
            subplot(length(preprocTypes), length(stimClasses), spNum); hold on;
            imagesc(strf); caxis([-ca ca]);
            title(sprintf('%s|%s: t=%0.2f, v=%0.2f', stimClass, preprocType, linTrainPerfs(k, m), linValidPerfs(k, m)));
            axis tight;
            
            figure(nlFig);
            subplot(length(preprocTypes), length(stimClasses), spNum); hold on;
            plot(x, fnval(dist.model.outputNL, x), 'r-', 'LineWidth', 2);
            plot(x, fnval(spline.model.outputNL, x), 'k-', 'LineWidth', 2);
            legend('Dist', 'Spline');
            title(sprintf('%s|%s: t=%0.2f -> (%0.2f vs %0.2f), v=%0.2f -> (%0.2f vs %0.2f)', stimClass, preprocType, ...
                  linTrainPerfs(k, m), distTrainPerfs(k, m), splineTrainPerfs(k, m), ...
                  linValidPerfs(k, m), distValidPerfs(k, m), splineValidPerfs(k, m)));
            axis([minX maxX 0 2]);
            
            figure(distFig);
            subplot(length(preprocTypes), length(stimClasses), spNum); hold on;
            errorbar(x, px_mean, px_std, 'k-');
            errorbar(x, pxspike_mean, pxspike_std, 'r-');
            errorbar(x, pxnospike_mean, pxnospike_std, 'b-');
            axis tight;
            
            figure(distFig2);
            subplot(length(preprocTypes), length(stimClasses), spNum); hold on;
            errorbar(x, pspikex_mean, pspikex_std, 'r-');
            errorbar(x, pnospikex_mean, pnospikex_std, 'b-');
            axis([min(x) max(x) 0 2]);
            
        end
        
    end
    
    distTrainDiffs = distTrainPerfs(:) - linTrainPerfs(:);
    distValidDiffs = distValidPerfs(:) - linValidPerfs(:);
    splineTrainDiffs = splineTrainPerfs(:) - linTrainPerfs(:);
    splineValidDiffs = splineValidPerfs(:) - linValidPerfs(:);
    
    figure(histFig);
    
    subplot(2, 2, 1); hold on;
    hist(distTrainDiffs);
    title('Dist: Train Diffs');
    subplot(2, 2, 2); hold on;
    hist(distValidDiffs);
    title('Dist: Valid Diffs');
    
    subplot(2, 2, 3); hold on;
    hist(splineTrainDiffs);
    title('Spline: Train Diffs');
    subplot(2, 2, 4); hold on;
    hist(splineValidDiffs);
    title('Spline: Valid Diffs');
        
end
    
    
function allData = get_all_unit_data(stimClasses, preprocTypes, outputDir)

    allData = struct;
    
    for k = 1:length(stimClasses)        
        stimClass = stimClasses{k};        
        allData.(stimClass) = struct;        
        for m = 1:length(preprocTypes)            
            preprocType = preprocTypes{m};
            
            pdata = struct;            
            responseFileBase = sprintf('directfit.%s.%s.h5', stimClass, preprocDesc);
            responseFileDists = sprintf('directfit.nl_dists.%s.%s.h5', stimClass, preprocDesc);
            responseFileSpline = sprintf('directfit.nl_spline.%s.%s.h5', stimClass, preprocDesc);

            mdata = get_model_data(fullfile(outputDir, responseFileBase));
            mdata_dists = get_model_data(fullfile(outputDir, responseFileDists));
            mdata_spline = get_model_data(fullfile(outputDir, responseFileSpline));
            
            pdata.linear = mdata;
            pdata.dists = mdata_dists;
            pdata.spline = mdata_spline;
            
            allData.(stimClass).(preprocType) = pdata;            
        end
    end
    
end