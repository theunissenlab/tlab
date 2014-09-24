function [stimFiles, respFiles] = pair_filenames_dcp(allStimsDir, respDir)


    stimLinkFiles = get_filenames(respDir, 'stim[0-9]*', 1);
    respFiles = get_filenames(respDir, 'spike[0-9]*', 1);
    
    stimFiles = cell(length(stimLinkFiles), 1);
    
    for k = 1:length(stimLinkFiles)       
        %read stim file and get path to .wav file
        stimLinkFile = stimLinkFiles{k};
        fid = fopen(stimLinkFile);
        wavFile = strtrim(fgetl(fid));
        fclose(fid);
        stimFiles{k} = fullfile(allStimsDir, wavFile);
    end
    