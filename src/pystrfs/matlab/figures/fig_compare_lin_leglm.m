function fig_compare_lin_leglm(cellName, preprocType, stimClass)

    descs = {'DirectFit', 'DirectFit+NL', 'Threshgrad', ...
             'Lin+Poisson', 'Lin+Binomial', ...             
             'LEGLM+Poisson'};

    fnames = {'directfit.%s.%s.%d.h5', 'directfit.nl_dists.%s.%s.%d.h5', ...
              'threshgrad.%s.%s.0.50.%d.h5', ...     
              'threshgrad.%s.%s.exponential.0.50.%d.h5', ...
              'threshgrad.%s.%s.logistic.0.50.%d.h5', ...              
              'threshgrad_glm.%s.%s.poisson.0.50.%d.h5'};

    data = get_modelcomp_cube(cellName, preprocType, stimClass, fnames);    
    infos = data.infos;
    strfs = data.strfs;
    biases = data.biases;
    
    figName = sprintf('%s | %s | %s', cellName, preprocType, stimClass);
    fprintf('\n%s\n', figName);
    fprintf('--------------------------\n');
    
    
    for k = 1:length(descs)
        perf = infos(k, :, 2) ./ infos(k, :, 1);
        fprintf('  %s:\n', descs{k});
        fprintf('    Performance: %0.3f +/- %0.3f\n', mean(perf), std(perf));
    end
       
    %% visualize strfs of TG, exp, log
    figure('Name', figName); hold on;
    dindx = 1:6;
    for k = 1:length(dindx)
        %indx = (k-1)*2 + 1;
        indx = k;
        di = dindx(k);
        
        strfMat = data.strfMats{di};
        
        if ~isempty(strfMat)
            subplot(3, 2, indx); hold on;
            strf = mean(strfMat, 3);
            bias = mean(biases(di, :));
            cval = max(abs(strf(:)));
            imagesc(strf); axis tight;
            caxis([-cval cval]);
            title(sprintf('%s: bias=%f', descs{di}, bias));
        end
        
    end
    