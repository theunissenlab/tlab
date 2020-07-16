function [stimFiles, respFiles] = pair_filenames(stimDir, respDir)
    
    if ~exist(stimDir, 'dir')
        error('Stimulus directory does not exist: %s\n', stimDir);
    end
    if ~exist(respDir, 'dir')
        error('Response directory does not exist: %s\n', respDir);
    end
    

    filePairs = struct;
    
    dlist = dir(respDir);
    pairNums = {};
    for k = 1:length(dlist)
        fname = dlist(k).name;
        respNum = regexp(fname, '[0-9]*', 'match');
        if ~isempty(respNum) && length(fname) > 2 && length(respNum) == 1
            respNum = ['s' num2str(respNum{1})];
            filePairs.(respNum) = struct;
            filePairs.(respNum).resp = fullfile(respDir, fname);
            pairNums{end+1} = respNum;
        end
    end
    
    dlist = dir(stimDir);
    for k = 1:length(dlist)       
        fname = dlist(k).name;        
        stimNum = regexp(fname, '[0-9]*', 'match');
        if ~isempty(stimNum) && length(fname) > 2 && length(stimNum) == 1
            stimNum = ['s' num2str(stimNum{1})];
            if isfield(filePairs, stimNum)
                filePairs.(stimNum).stim = fullfile(stimDir, fname);
            end
        end
    end
    
    pairCount = length(pairNums);
    stimFiles = cell(pairCount, 1);
    respFiles = cell(pairCount, 1);
    
    for k = 1:length(pairNums)        
        stimFiles{k} = filePairs.(pairNums{k}).stim;
        respFiles{k} = filePairs.(pairNums{k}).resp;
    end
    