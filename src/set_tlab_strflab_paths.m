function set_tlab_strflab_paths(tlabTrunk, strflabTrunk)

    tlabDirs = {'file', 'spike', 'sound', 'STRF', 'direct_fit', 'validation', 'strflab', 'strflab/auxFunc', 'strflab/layer_3optim', 'strflab/models/lnp', 'strflab/models/lnpsr', 'strflab/models/lnpsr2', '../AuditoryToolbox'};
    strflabDirs = {'.', 'auxFunc', 'auxFunc/direct_fit', 'layer_1', 'layer_2resamp', 'layer_3optim', 'models/lin', 'preprocessing', 'preprocessing/sound', 'validation'};
    
    %add strflab directories first
    for k = 1:length(strflabDirs)       
        dname = fullfile(strflabTrunk, strflabDirs{k});
        addpath(dname);        
    end
    
    %add tlab directories so duplicate files are read from here
    for k = 1:length(tlabDirs)       
        dname = fullfile(tlabTrunk, tlabDirs{k});
        addpath(dname);        
    end
       