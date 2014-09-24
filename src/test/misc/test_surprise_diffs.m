function test_surprise_diffs()

    oldCell = 'yy1617_4_A';
    oldDataDir = fullfile('/auto/fdata/mschachter/data', oldCell, 'conspecific');
    oldDataFile = fullfile(oldDataDir, 'preproc', 'surprise.stft.default.mat');
    md5s = get_stim_md5s(oldDataDir);
    
    preprocFile = '/auto/k6/mschachter/pystrfs/preproc/surprise.dfw_3.dtw_3.dg_4.Con.h5';
    [newStims, newMd5Map] = get_stim_data(preprocFile);
    
    preprocFile = '/auto/k6/mschachter/pystrfs/preproc/stft.nstd_6.fband_125.h5';
    [stftStims, stftMd5Map] = get_stim_data(preprocFile);
    
    vars = load(oldDataFile);
    oldStimLouder = vars.surpriseStimLouder;
    oldStimQuieter = vars.surpriseStimQuieter;
    oldStimFull = [oldStimLouder; oldStimQuieter];
    groupIndex = vars.groupIndex;
    
    stimPairs = cell(length(md5s), 3);
    
    for k = 1:length(md5s)        
        nsIndx = newMd5Map.(['m' md5s{k}]);
        newStim = newStims{nsIndx}.spec;   
        
        stIndx = stftMd5Map.(['m' md5s{k}]);
        stftStim = stftStims{stIndx}.spec;
        
        oldStim = oldStimFull(:, find(groupIndex == k));
        
        stimPairs{k, 1} = oldStim;
        stimPairs{k, 2} = newStim;        
        stimPairs{k, 3} = stftStim;
    end
    
    
    for k = 1:length(md5s)
        
        oldStim = stimPairs{k, 1};
        alen = round(size(oldStim, 1) / 2);
        oldLouder = oldStim(1:alen, :);
        oldQuieter = oldStim(alen+1:end, :);
                
        newStim = stimPairs{k, 2};
        newLouder = newStim(1:alen, :);
        newQuieter = newStim(alen+1:end, :);
        
        
        figure; hold on;
                
        h = subplot(3, 2, 1); hold on;
        p = get(h, 'pos');
        p(3) = p(3) + 0.13;
        p(4) = p(4) + 0.03;
        %set(h, 'pos', p);
        imagesc(stimPairs{k, 3}); axis tight; colorbar;
        title('STFT');
        
        h = subplot(3, 2, 2); hold on;
        p = get(h, 'pos');
        p(3) = p(3) + 0.13;
        p(4) = p(4) + 0.03;
        %set(h, 'pos', p);
        imagesc(stimPairs{k, 3}); axis tight; colorbar;
        title('STFT');
        
        h = subplot(3, 2, 3); hold on;
        p = get(h, 'pos');
        p(3) = p(3) + 0.13;
        p(4) = p(4) + 0.03;
        %set(h, 'pos', p);
        imagesc(oldLouder); axis tight; colorbar;
        title('Old Louder');
        
        h = subplot(3, 2, 5); hold on;
        p = get(h, 'pos');
        p(3) = p(3) + 0.13;
        p(4) = p(4) + 0.03;
        %set(h, 'pos', p);
        imagesc(newLouder); axis tight; colorbar;
        title('New Louder');
        
        h = subplot(3, 2, 4); hold on;
        p = get(h, 'pos');
        p(3) = p(3) + 0.13;
        p(4) = p(4) + 0.03;
        %set(h, 'pos', p);
        imagesc(oldQuieter); axis tight; colorbar;
        title('Old Quieter');
        
        h = subplot(3, 2, 6); hold on;
        p = get(h, 'pos');
        p(3) = p(3) + 0.13;
        p(4) = p(4) + 0.03;
        %set(h, 'pos', p);
        imagesc(newQuieter); axis tight; colorbar;
        title('New Quieter');
        
        
    end
    
end



function md5s = get_stim_md5s(dataDir, stimDir)

    stimLinkFiles = get_filenames(dataDir, 'stim[0-9]*', 1);
    
    if iscell(stimLinkFiles)

        nFiles = length(stimLinkFiles);        
        md5s = cell(nFiles, 1);
        
        for k = 1:nFiles
            
            %read stim file and get path to .wav file
            stimLinkFile = stimLinkFiles{k};
            fid = fopen(stimLinkFile);
            wavFile = strtrim(fgetl(fid));
            fclose(fid);
            
            [rootDir, md5, ext] = fileparts(wavFile);
            
            md5s{k} = md5;
        end
        
    end
end