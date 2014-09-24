function preprocess_spectrogram(stimFiles, md5s, outputFile, tfParams)
    
%               .fband: spacing between frequency bands of spectrogram (125)
%               .nstd: # of std deviations that define width of Gaussian
%                   window (6)
%               .low_freq: lowest frequency in Hz (250)
%               .high_freq: highest frequency in Hz (8000)
%               .log: take base 10 log of spectrogram (1)
%               .refpow: reference power when taking log spectrogram
%               .dbnoise: noise floor in dB, everything below this is
%                   ignored ( when .log=1)

    %% preprocess the sound to produce spectrograms
    preprocStimParams = struct;      
    preprocStimParams.tfType = 'stft';
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
        h5.set_attr(fid, groupName, 'fband', preprocStimParams.tfParams.fband);
        h5.set_attr(fid, groupName, 'nstd', preprocStimParams.tfParams.nstd);
        h5.set_attr(fid, groupName, 'low_freq', preprocStimParams.tfParams.low_freq);
        h5.set_attr(fid, groupName, 'high_freq', preprocStimParams.tfParams.high_freq);
        h5.set_attr(fid, groupName, 'log', preprocStimParams.tfParams.log);
        h5.set_attr(fid, groupName, 'refpow', preprocStimParams.tfParams.refpow);
        h5.set_attr(fid, groupName, 'dbnoise', preprocStimParams.tfParams.dbnoise);
        h5.set_attr(fid, groupName, 'sample_rate', stimInfo.sampleRate);    
        h5.set_attr(fid, groupName, 'max_power', maxPow);    

    end
    
    H5F.close(fid);  