function lnl_comp_onecell_display(onecell_data, ptypes, stypes)

    if nargin < 2
        ptypes = onecell_data.preprocTypes;
    end
    if nargin < 3
        stypes = onecell_data.stimTypes;
    end

    ncols = length(stypes);
    npreprocs = length(ptypes);
        
    %% figure 1: STRFS
    figure('name', sprintf('STRFs: %s', onecell_data.cellName));
    for k = 1:length(ptypes)
        for m = 1:length(stypes)
            
            ptype = ptypes{k};
            stype = stypes{m};       
            odata = onecell_data.(ptype).(stype);
            indx = m + (k-1)*ncols;
            subplot(npreprocs, ncols, indx); hold on;
            
            strf = odata.nl.strf;
            b = odata.nl.bias;
            mval = max(abs(strf(:)));
            linPerf = odata.lin.perf;
            nlPerf = odata.nl.perf;            
            
            imagesc(strf); axis tight;
            caxis([-mval mval]); colorbar;
            title(sprintf('%s | %s | b=%0.4f | lin=%0.2f, nl=%0.2f', ...
                  ptypes{k}, stypes{m}, b, linPerf, nlPerf));
        end
        
    end

    %% figure 3: Diffs between STRFs
    figure('name', sprintf('STRF Differences: %s', onecell_data.cellName));
    for k = 1:length(ptypes)
        for m = 1:length(stypes)
            
            ptype = ptypes{k};
            stype = stypes{m};       
            odata = onecell_data.(ptype).(stype);
            indx = m + (k-1)*ncols;
            subplot(npreprocs, ncols, indx); hold on;
            
            origStrf = odata.nl.origStrf;
            nlStrf = odata.nl.strf;
            origBias = odata.nl.origBias;
            nlBias = odata.nl.bias;
            bestIter = odata.nl.bestIter;
            strfDiff = nlStrf - origStrf;
            mval = max(abs(strfDiff(:)));
                        
            imagesc(strfDiff); axis tight;
            caxis([-mval mval]); colorbar;
            title(sprintf('%s | %s | bias: %0.4f->%0.4f | bestiter=%d', ...
                  ptypes{k}, stypes{m}, origBias, nlBias, bestIter));
        end
        
    end
    
    %% figure 2: Output NLs
    figure('name', sprintf('Fit Output NL: %s', onecell_data.cellName));
    for k = 1:length(ptypes)
        for m = 1:length(stypes)
            
            ptype = ptypes{k};
            stype = stypes{m};       
            odata = onecell_data.(ptype).(stype);
            indx = m + (k-1)*ncols;
            subplot(npreprocs, ncols, indx); hold on;
            
            maxResp = odata.nl.maxLinResp;
            minResp = odata.nl.minLinResp;
            nlInc = (maxResp - minResp) / 100;
            x = minResp:nlInc:maxResp;

            b = odata.nl.b;
            c = odata.nl.c;
            linPerf = odata.lin.perf;
            nlPerf = odata.nl.perf;            
            
            y = gompertz(1, b, c, x);
            plot(x, y, 'k-', 'LineWidth', 2);
            errorbar(odata.nl.binned.x, odata.nl.binned.y, odata.nl.binned.stds);
            axis tight;
            title(sprintf('%s | %s | b=%0.1f | c=%0.1f | lin=%0.2f, nl=%0.2f', ...
                  ptypes{k}, stypes{m}, b, c, linPerf, nlPerf));            
        end        
    end
    
    
    %% figure 4: Response Comparisons
    figure('name', sprintf('Responses: %s', onecell_data.cellName));
    for k = 1:length(ptypes)
        for m = 1:length(stypes)
            
            ptype = ptypes{k};
            stype = stypes{m};       
            odata = onecell_data.(ptype).(stype);
            indx = m + (k-1)*ncols;
            subplot(npreprocs, ncols, indx); hold on;
            
            linPerf = odata.lin.perf;
            nlPerf = odata.nl.perf;            
            
            plot(odata.realResp, 'k-', 'LineWidth', 2);
            plot(odata.modelRespLin, 'r-');
            plot(odata.modelRespNL, 'b-');
            axis tight;
            
            title(sprintf('%s | %s | lin=%0.2f, nl=%0.2f', ...
                  ptypes{k}, stypes{m}, linPerf, nlPerf));
        end
        
    end    
    
    
    %% figure 5: Response Differences
    figure('name', sprintf('Response Differences: %s', onecell_data.cellName));
    for k = 1:length(ptypes)
        for m = 1:length(stypes)
            
            ptype = ptypes{k};
            stype = stypes{m};       
            odata = onecell_data.(ptype).(stype);
            indx = m + (k-1)*ncols;
            subplot(npreprocs, ncols, indx); hold on;
            
            linPerf = odata.lin.perf;
            nlPerf = odata.nl.perf;            
            
            
            linDiff = rv(odata.modelRespLin) - rv(odata.realResp);
            nlDiff = rv(odata.modelRespNL) - rv(odata.realResp);
            
            plot(linDiff, 'r-');
            plot(nlDiff, 'b-');
            axis tight;
            
            title(sprintf('%s | %s | lin=%0.2f, nl=%0.2f', ...
                  ptypes{k}, stypes{m}, linPerf, nlPerf));
            
        end
    end
    
    
    %% figure 6: Modspecs of STRFs
    
    