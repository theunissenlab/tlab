function stim_sample_rates()
    
    stimDir = '/auto/k6/tlab/crcns/all_stims';
    
    fileNames = get_filenames(stimDir, '.wav', 1);
    
    if ~iscell(fileNames)
        fprintf('No filenames found in %s...\n', stimDir);
        return;
    end
    
    fout = fopen('/auto/k6/tlab/crcns/stim_data.csv', 'w');
    for k = 1:length(fileNames)
        
        fname = fileNames{k};               
        [pstr, name, ext] = fileparts(fname);
        [adata, sampleRate, depth] = wavread(fname);        
        alen = length(adata);
        fprintf(fout, '%s,%f,%d,%d\n', [name ext], sampleRate, depth, alen);
        clear adata;
        
    end
    fclose(fout);