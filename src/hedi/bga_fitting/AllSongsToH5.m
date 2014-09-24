function AllSongsToH5(dirname, outputFile)
    %% takes a directory filled with motogram .mat files and writes them to an hdf5 file
    
    fnames = dir(dirname);
    
    h5 = h5utils();
    f = h5.create(outputFile);
    
    motogramGroup = '/motograms';
    h5.create_group(f, motogramGroup);
    
    for k = 1:length(fnames)
        
        fname = fnames(k).name;
        if strcmp(fname(1), '.')
            continue;
        end
        fullFileName = fullfile(dirname, fname);
        
        vars = load(fullFileName);
        SongToH5(vars.song, f, motogramGroup);
    end
    
    h5.close(f);