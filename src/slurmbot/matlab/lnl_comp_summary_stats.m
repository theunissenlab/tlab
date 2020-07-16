function stats = lnl_comp_summary_stats(alldata)

    indx.mld = 1:9;
    indx.ov = 10:19;
    indx.fldl = 20:29;

    regions = {'mld', 'ov', 'fldl'};
    preprocTypes = {'stft', 'lyons', 'surprise'};
    stimTypes = {'conspecific', 'flatrip'};
    outputNLTypes = {'lin', 'nl'};
    
    nstim = length(stimTypes);
    npre = length(preprocTypes);
    nout = length(outputNLTypes);
    n = length(alldata);
    
    stats = struct;
    stats.perf = zeros(n, npre, nstim, nout);
    stats.infoBound = zeros(n, npre, nstim);
    
    
    %% put data into hypercube
    for k = 1:length(alldata)    
        onecell_data = alldata{k};        
        for i = 1:length(preprocTypes)
            for j = 1:length(stimTypes)
                odata = onecell_data.(preprocTypes{i}).(stimTypes{j});
                stats.infoBound(k, i, j) = odata.infoBound;
                for m = 1:length(outputNLTypes)
                    stats.perf(k, i, j, m) = odata.(outputNLTypes{m}).perf;                    
                end
            end
        end        
    end
    
    
    stimType = 1;    
    %% analyze performance histograms for sub-classes
    for k = 1:length(regions)
        for i = 1:length(preprocTypes)
            
            regIndx = indx.(regions{k});                

            linPerfs = stats.perf(regIndx, i, stimType, 1);
            nlPerfs = stats.perf(regIndx, i, stimType, 2);
            perfDiffs = nlPerfs - linPerfs;

            [h1, p1, ci1] = ttest(perfDiffs);

            pdMean = mean(perfDiffs);
            pdStd = std(perfDiffs);
            ttl = sprintf('%s %s %s  |  %0.3f +/- %0.4f  | nonzero-mean=%d; p=%0.2f [%0.3f, %0.3f]', ...
                          regions{k}, preprocTypes{i}, stimTypes{stimType}, ...
                          pdMean, pdStd, ...
                          h1, p1, ci1(1), ci1(2));
            fprintf(ttl);
            [linPerfs nlPerfs perfDiffs]

            %{
            figure('name', 'Perf Differences'); hold on;
            hist(perfDiffs);                
            title(ttl);
            %}

        end        
    end
    
    %{
    linPerfsCon = stats.perf(:, 1, 1, 1);
    nlPerfsCon = stats.perf(:, 1, 1, 2);
    perfDiffsCon = nlPerfsCon - linPerfsCon;
    
    figure; hold on;
    hist(perfDiffsCon, 15);
    title('Spectrogram + Conspecific');
    
    
    linPerfsCon = stats.perf(:, 3, 1, 1);
    nlPerfsCon = stats.perf(:, 3, 1, 2);
    perfDiffsCon = nlPerfsCon - linPerfsCon;
    
    figure; hold on;
    hist(perfDiffsCon, 15);
    title('Surprise + Conspecific');
    %}
    
    
    for m = 1:length(regions)
       
        reg = regions{m};
       
        ncells = 6;
        figure('name', reg); hold on;
        for k = 1:ncells       
            cindx = indx.(reg);
            c = cindx(k);
            odata = alldata{c}.stft.conspecific;
            
            spindx = (k-1)*2 + 1;
            
            strf = odata.nl.strf;
            subplot(ncells, 2, spindx); hold on;
            mstrf = max(abs(strf(:)));
            imagesc(strf); caxis([-mstrf mstrf]); axis tight;
            
            subplot(ncells, 2, spindx+1); hold on;
            errorbar(odata.nl.binned.x, odata.nl.binned.y, odata.nl.binned.stds);
            axis tight;            
            
        end
        
    end
        
        
        
        
    end
    
    
    
    
    



    