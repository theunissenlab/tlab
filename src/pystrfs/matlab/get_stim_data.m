function [stimData, md5Map] = get_stim_data(preprocFile)

    h5 = h5utils();
    
    pfid = h5.open(preprocFile);    
    
    stimMd5s = h5.get_subgroups(pfid, '/');    
    stimData = cell(length(stimMd5s), 1);    
    
    md5Map = struct;
    
    for k = 1:length(stimMd5s)
    
        stim = struct;
        
        %% get spectrogram
        specGroup = sprintf('/%s', stimMd5s{k});
        stimDs = sprintf('%s/spectrogram', specGroup);
        
        stim.lowFreq = h5.get_attr(pfid, specGroup, 'low_freq');
        stim.highFreq = h5.get_attr(pfid, specGroup, 'high_freq');
        stim.sampleRate = h5.get_attr(pfid, specGroup, 'sample_rate');
        
        stim.spec = h5.get_ds(pfid, stimDs);        
        stim.length = size(stimData{k}, 2) / stim.sampleRate;
        stim.md5 = stimMd5s{k};
        
        stimData{k} = stim;
        
        md5Map.(['m' stimMd5s{k}]) = k;        
    end
    