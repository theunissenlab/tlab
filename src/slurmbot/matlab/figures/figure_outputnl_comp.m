function figure_outputnl_comp(cellNames)
    
    cinfo = get_cell_info();
    
    for k = 1:length(cellNames)
        
        cname = cellNames{k};        
        stypes = cinfo.(cname).stimTypes;        
        nstypes = length(stypes);
        
        figure('Name', sprintf('%s | %s', cname, cinfo.(cname).region)); hold on;
        
        maxY = -1;
        xBnds = zeros(nstypes*2, 2);
        for m = 1:nstypes
        
            stimType = stypes{m};
            outputDir = fullfile(cinfo.(cname).rootDir, stimType, 'output');
                       
            subplot(nstypes, 2, m*2 - 1); hold on;
            stftName = sprintf('strflab.tfType_stft.%s.mat', stimType);
            stftFullPath = fullfile(outputDir, stftName);            
            [xvals, yvals, stds, perfData, gps] = display_outputnl_info(stftFullPath);
            gvals = gompertz(gps(1), gps(2), gps(3), xvals);
            
            xBnds(m*2-1, 1) = min(xvals);
            xBnds(m*2-1, 2) = max(xvals);
            maxY = max([maxY cv(yvals) cv(gvals)]);
            errorbar(xvals, yvals, stds);            
            plot(xvals, gvals, 'r-');
            
            title(sprintf('stft | %s | %0.3f | %0.3f', stimType, perfData.train.perf, perfData.valid.perf));
            
            subplot(nstypes, 2, m*2); hold on;
            lyonsFullPath = get_best_lyons(outputDir);
            [xvals, yvals, stds, perfData, gps] = display_outputnl_info(lyonsFullPath);
            gvals = gompertz(gps(1), gps(2), gps(3), xvals);
            
            errorbar(xvals, yvals, stds);
            plot(xvals, gvals, 'r-');
            title(sprintf('lyons | %s | %0.3f | %0.3f', stimType, perfData.train.perf, perfData.valid.perf));
            maxY = max([maxY cv(yvals) cv(gvals)]);
            xBnds(m*2, 1) = min(xvals);
            xBnds(m*2, 2) = max(xvals);
        end        
        
        for m = 1:(nstypes*2)
            subplot(nstypes, 2, m); hold on;
            axis([xBnds(m, 1) xBnds(m, 2) 0 maxY]);
        end
    end
end
    
function outputFile = get_best_lyons(outputDir)   
    bfile = fullfile(outputDir, 'best.lyons.txt');    
    fid = fopen(bfile);
    a = textscan(fid, '%s');
    outputFile = a{1}{1};
    outputFile = strrep(outputFile, '.txt', '.mat');
end
