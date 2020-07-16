function unit = read_unit_file(unitPath)

    h5 = h5utils();

    unit = struct;
    
    fid = h5.open(unitPath, 'r');
    
    unit.description = h5.get_attr(fid, '/', 'description');
    unit.protocol_id = h5.get_attr(fid, '/', 'protocol_id');
    unit.unit_id = h5.get_attr(fid, '/', 'unit_id');
    
    crList = h5.get_subgroups(fid, '/class_responses');
    unit.classes = crList;
    
    for k = 1:length(crList)

        crPath = sprintf('/class_responses/%s', crList{k});
        numStims = h5.get_attr(fid, crPath, 'num_stims');
        
        resps = {};
        
        for m = 1:int8(numStims)
           
            respPath = sprintf('%s/%d', crPath, m);
            
            numTrials = h5.get_attr(fid, respPath, 'num_trials');
            protocolStimId = h5.get_attr(fid, respPath, 'protocol_stim_id');
            trialDuration = h5.get_attr(fid, respPath, 'trial_duration');
            stimMd5 = h5.get_attr(fid, respPath, 'stim_md5');
            
            spikeTrials = {};
            for p = 1:int8(numTrials)
                stPath = sprintf('%s/%d', respPath, p);
                stimes = h5.get_ds(fid, stPath);
                spikeTrials{end+1} = stimes;                
            end            
            
            resp = struct;
            resp.numTrials = numTrials;
            resp.protocolStimId = protocolStimId;
            resp.spikeTrials = spikeTrials;
            resp.trialDuration = trialDuration;
            resp.stimMd5 = stimMd5;
            
            resps{end+1} = resp;
        end
        
        unit.class_responses.(crList{k}) = resps;
        
    end
    
    h5.close(fid);
end
