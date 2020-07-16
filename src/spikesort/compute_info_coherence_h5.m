function compute_info_coherence_h5(unitFile, outputFile)

    unit = read_unit_h5file(unitFile);
    
    sampleRate = 1000;
    freqCutoff = 500;
    
    h5 = h5utils();
    fid = h5.create(outputFile);
        
    for k = 1:length(unit.classes)
    
        className = unit.classes{k};
        resps = unit.class_responses.(className);
        
        %% combine all responses within class into cell array
        allSpikeTrials = cell(length(resps), 1);
        stimLengths = zeros(length(resps), 1);
        for m = 1:length(resps)
            rtrials = resps{m}.trials;
            ntrials = length(rtrials);
            stimes = cell(ntrials, 1);
            for n = 1:ntrials
                stimes{n} = rtrials{n}.spikeTimes;                
            end
            allSpikeTrials{m} = stimes;
            stimLengths(m) = resps{m}.stim_duration;
        end

        preprocRespParams = struct;         %create preprocessing struct
        preprocRespParams.units = 's';      %files have units of seconds
        preprocRespParams.binSize = 0.001;  %1ms bin size, same sample rate as spectrograms
        preprocRespParams.stimLengths = stimLengths; %needed to compute correct PSTH lengths
        
        [wholeResponse, groupIndex, respInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);


        %% split original spike trials into two PSTHs for purposes of validation
        preprocRespParams.split = 1;
        [wholeSplitResponse, respSplitInfo, preprocRespParams] = preprocSpikeResponses(allSpikeTrials, preprocRespParams);

        avgNumTrials = mean(respInfo.numTrials);

        infoWindowSize = 0.200;

        cStruct = compute_coherence_bound(wholeSplitResponse(1, :), wholeSplitResponse(2, :), avgNumTrials, sampleRate, freqCutoff, infoWindowSize);        
                    
        meanSpikeRate = mean(wholeResponse) * sampleRate;        
        
        pathPrefix = sprintf('/class_info/%s', className);
        h5.set_ds(fid, pathPrefix, 'mean_rate', meanSpikeRate);
        h5.set_ds(fid, [pathPrefix '/coherence'], 'cutoff_frequency', cStruct.freqCutoff);
        h5.set_ds(fid, [pathPrefix '/coherence'], 'frequencies', cStruct.f);
        h5.set_ds(fid, [pathPrefix '/coherence'], 'upper', cStruct.cUpper);
        h5.set_ds(fid, [pathPrefix '/coherence'], 'lower', cStruct.cLower);
        h5.set_ds(fid, [pathPrefix '/coherence'], 'mean', cStruct.c);
        
        h5.set_ds(fid, [pathPrefix '/coherence/info'], 'upper', cStruct.infoUpper);
        h5.set_ds(fid, [pathPrefix '/coherence/info'], 'lower', cStruct.infoLower);
        h5.set_ds(fid, [pathPrefix '/coherence/info'], 'mean', cStruct.info);        
        
    end
    
    h5.close(fid);
        
        
end