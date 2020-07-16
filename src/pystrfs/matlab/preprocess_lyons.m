function preprocess_lyons(stimFiles, md5s, outputFile, tfParams)

%           'lyons' params (Requires AuditoryToolbox)               
%               .low_freq: lowest frequency in Hz (250)
%               .high_freq: highest frequency in Hz (8000)
%               .earQ: quality factor of each filter (8)
%               .agc: use adaptive gain control (1)
%               .differ: use differential gain control (1)
%               .tau: time constant of gain control (3)
%               .step: 1/step is approximately the number of
%               filters per bandwidth 


    %% preprocess the sound to produce Lyons cochleagrams
    preprocStimParams = struct;      
    preprocStimParams.tfType = 'lyons';
    preprocStimParams.tfParams = tfParams;
    preprocStimParams.cache = 0;
    
    if exist(outputFile, 'file')
        delete(outputFile);
    end
    
    h5 = h5utils();
    fid = h5.create_file(outputFile);
    
    for k = 1:length(stimFiles)    
        fprintf('Preprocessing %s...\n', stimFiles{k});
        [spec, groupIndex, stimInfo, preprocStimParams] = preprocSound({stimFiles{k}}, preprocStimParams);
        
        maxPow = max(spec(:));
        
        groupName = sprintf('/%s', md5s{k});
        h5.set_ds(fid, groupName, 'spectrogram', spec');
	specSize = size(spec)

        h5.set_attr(fid, groupName, 'length', stimInfo.stimLengths(1));    
        h5.set_attr(fid, groupName, 'earQ', preprocStimParams.tfParams.earQ);
        h5.set_attr(fid, groupName, 'agc', preprocStimParams.tfParams.agc);
        h5.set_attr(fid, groupName, 'differ', preprocStimParams.tfParams.differ);
        h5.set_attr(fid, groupName, 'low_freq', preprocStimParams.tfParams.low_freq);
        h5.set_attr(fid, groupName, 'high_freq', preprocStimParams.tfParams.high_freq);
        h5.set_attr(fid, groupName, 'log', preprocStimParams.tfParams.log);
        h5.set_attr(fid, groupName, 'tau', preprocStimParams.tfParams.tau);
        h5.set_attr(fid, groupName, 'step', preprocStimParams.tfParams.step);
        h5.set_attr(fid, groupName, 'sample_rate', stimInfo.sampleRate);    
        h5.set_attr(fid, groupName, 'max_power', maxPow);    

    end
    
    H5F.close(fid);  