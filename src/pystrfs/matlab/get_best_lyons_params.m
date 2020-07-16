function [earQ, step] = get_best_lyons_params(cellName, stimClass)

    outputDir = fullfile('/auto/k6/mschachter/pystrfs/units', cellName, 'output');
    
    bestFileName = sprintf('lyons.bestparams.%s.csv', stimClass);
    bestFile = fullfile(outputDir, bestFileName);
    
    fid = fopen(bestFile, 'r');
    ln = fgetl(fid);
    while ischar(ln)
        
        toks = regexp(ln, '=', 'split');
        if length(toks) > 0           
            vname = toks{1};
            vval = toks{2};
            if strcmp(vname, 'earQ')
                earQ = str2double(vval);
            elseif strcmp(vname, 'step')
                step = str2double(vval);
            end
        end
        ln = fgetl(fid);
    end
    