function compute_neurogram_coherence_for_all(max_id)

    root_dir = '/auto/k6/mschachter/netman/neurograms';
        
    for k = 1:max_id
        
        ng_file = fullfile(root_dir, sprintf('neurogram_response_%d.h5', k));
        strf_file = fullfile(root_dir, sprintf('strf_response_%d.h5', k));
        
        if exist(ng_file, 'file') && exist(strf_file, 'file')
            compute_response_info_coherence(ng_file);
            compute_response_info_coherence(strf_file);            
        end        
    end
    
    
    