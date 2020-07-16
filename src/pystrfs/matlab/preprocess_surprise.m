function preprocess_surprise(md5s, preprocessFile, params, transforms, outputFile)
    
    surpriseParams = struct;
    surpriseParams.domainFrequencyWidth = params.domainFrequencyWidth;
    surpriseParams.domainTimeWidth = params.domainTimeWidth;
    surpriseParams.domainGap = params.domainGap;
    surpriseParams.cache = 0;
    
    if exist(outputFile, 'file')
        delete(outputFile);
    end

    %% load preprocessed stimuli and concatenate them
    h5 = h5utils();
    fid = h5.open(preprocessFile);
    preprocStims = cell(length(md5s), 1);
    stimInfos = cell(length(md5s), 1);
    numChannels = -1;
    totalLen = 0;
    groupIndex = [];
    for k = 1:length(preprocStims)
      
        groupName = sprintf('/%s', md5s{k});
        dsPath = sprintf('%s/spectrogram', groupName);
        preprocStims{k} = h5.get_ds(fid, dsPath);             
        if numChannels == -1
            numChannels = size(preprocStims{k}, 1);
        end
        tlen = size(preprocStims{k}, 2);
        totalLen = totalLen + tlen;
        groupIndex = [groupIndex ones(1, tlen)*k];
        
        sinfo = struct;
        sinfo.length = h5.get_attr(fid, groupName, 'length');
        sinfo.low_freq = h5.get_attr(fid, groupName, 'low_freq');
        sinfo.high_freq = h5.get_attr(fid, groupName, 'high_freq');
        sinfo.sample_rate =  h5.get_attr(fid, groupName, 'sample_rate');
        stimInfos{k} = sinfo;
        
    end
    h5.close(fid);

    wholeStim = zeros(totalLen, numChannels);
    for k = 1:length(preprocStims)
        tindx = find(groupIndex==k);
        wholeStim(tindx, :) = preprocStims{k}';
    end
    wholeStimSize = size(wholeStim)
    groupIndexSize = size(groupIndex)
   
    clear preprocStims;
    
    for k = 1:length(transforms)       
        fprintf('Performing %s transform...\n', transforms{k});
        tfunc = sprintf('transform_%s(wholeStim)', transforms{k});
        wholeStim = eval(tfunc);        
    end
    
    
    fprintf('Preprocessing surprise for %d stimuli...\n', length(md5s));
    %% preprocess surprise
    [surpriseStimLouder, surpriseStimQuieter, groupIndex, stimInfo, surpriseParams] = preprocSoundSurprise(wholeStim, groupIndex, surpriseParams);    
    
    clear wholeStim;
    
    h5 = h5utils();
    fid = h5.create(outputFile);

    for k = 1:length(md5s)
       
        indx = find(groupIndex == k);                
        spec = [surpriseStimLouder(:, indx); surpriseStimQuieter(:, indx)];        
        
        groupName = sprintf('/%s', md5s{k});
        h5.set_ds(fid, groupName, 'spectrogram', spec);

        sinfo = stimInfos{k};
        h5.set_attr(fid, groupName, 'length', sinfo.length);            
        h5.set_attr(fid, groupName, 'low_freq', sinfo.low_freq);
        h5.set_attr(fid, groupName, 'high_freq', sinfo.high_freq);
        h5.set_attr(fid, groupName, 'sample_rate', sinfo.sample_rate);            
        
    end

    h5.close(fid);  
    