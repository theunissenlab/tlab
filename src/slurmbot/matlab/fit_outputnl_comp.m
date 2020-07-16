function fit_outputnl_comp(cellNames)
    
    cinfo = get_cell_info();
    
    for k = 1:length(cellNames)
        
        cname = cellNames{k};        
        stypes = cinfo.(cname).stimTypes;        
        nstypes = length(stypes);
        
        for m = 1:nstypes
        
            stimType = stypes{m};
            outputDir = fullfile(cinfo.(cname).rootDir, stimType, 'output');
                       
            stftName = sprintf('strflab.tfType_stft.%s.mat', stimType);
            stftFullPath = fullfile(outputDir, stftName);            
            [xvals, yvals, stds, perfData, gps] = display_outputnl_info(stftFullPath);
            
            fname = 'outputNL.gompertz.tfType_stft.txt';
            ofile = fullfile(outputDir, fname);
            f = fopen(ofile, 'w');
            fprintf(f, '%0.6f,%0.6f,%0.6f\n', gps(1), gps(2), gps(3));
            fclose(f);
            
            lyonsFullPath = get_best_lyons(outputDir);
            [xvals, yvals, stds, perfData, gps] = display_outputnl_info(lyonsFullPath);
            
            fname = 'outputNL.gompertz.tfType_lyons.txt';
            ofile = fullfile(outputDir, fname);
            f = fopen(ofile, 'w');
            fprintf(f, '%0.6f,%0.6f,%0.6f\n', gps(1), gps(2), gps(3));
            fclose(f);            
            
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
