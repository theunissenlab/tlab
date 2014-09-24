function batch_dataset_coherence(rootdir, stimdir)
     
    datasets = find_datasets(rootdir, stimdir, 'conspecific|zfsongs');

    for k = 1:length(datasets)
        tic;
        ds = datasets{k};
        
        dataDir = ds.dirname;
        stimFiles = ds.srPairs.stimFiles;
        respFiles = ds.srPairs.respFiles;
        
        compute_coherence_dataset(dataDir, stimFiles, respFiles, 0);
        
        etime = toc;
        fprintf('Took %f seconds, %d to go...\n', etime, (length(datasets)-k));        
    end
    