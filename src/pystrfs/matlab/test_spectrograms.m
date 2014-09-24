function test_spectrograms()


    dataDir = '/auto/fdata/mschachter/data';
    allStimsDir = fullfile(dataDir, 'all_stims');
    cellName = 'yy1617_4_A';
    cellDir = fullfile(dataDir, cellName);

    stimDir = fullfile(cellDir, 'conspecific');
        
    srPairs = get_sr_files(allStimsDir, stimDir);
    stimFiles = cell(length(srPairs.stimFiles), 1);
    md5s = cell(length(srPairs.stimFiles), 1);
    for k = 1:length(md5s)      
      [path,name,ext] = fileparts(srPairs.stimFiles{k});
      md5s{k} = name;
      stimFiles{k} = srPairs.stimFiles{k};
    end

    %% get data from old preprocessing
    preprocFiles = cell(length(md5s), 1);
    preprocSpecs = cell(length(md5s), 1);
    for k = 1:length(preprocFiles)
      pname = sprintf('stim.preproc-stft.default.%d.mat', k);
      pfile = fullfile(stimDir, 'preproc', pname);
      preprocFiles{k} = pfile;
      pvars = load(pfile);
      spec = pvars.stim.tfrep.spec;
      preprocSpecs{k} = spec;
    end

    h5 = h5utils();
    fid = H5F.open('/auto/k6/mschachter/pystrfs/preproc/stft.nstd_6.fband_125.h5',...
		   'H5F_ACC_RDONLY','H5P_DEFAULT');
    h5Specs = cell(length(md5s), 1);
    %% get data from h5 file
    for k = 1:length(md5s)

      groupPath = sprintf('/%s/spectrogram', md5s{k});
      spec = h5.get_ds(fid, groupPath);
      logSpec = max(0, 20*log10(spec/6)+80);
      h5Specs{k} = logSpec;
    end

    H5F.close(fid);
    
    ofname = sprintf('/auto/fhome/mschachter/spec_comp_%s.mat', ...
		     cellName);
    save(ofname, 'preprocSpecs', 'h5Specs');

    

    
function srPairs = get_sr_files(allStimsDir, dataDir)

    srPairs = -1;
    stimLinkFiles = get_filenames(dataDir, 'stim[0-9]*', 1);
    respFiles = get_filenames(dataDir, 'spike[0-9]*', 1);
    
    if iscell(stimLinkFiles) && iscell(respFiles) && (length(stimLinkFiles) ...
						      == ...
						      length(respFiles))
       
        srPairs = struct;
        nFiles = length(stimLinkFiles);        
        stimFiles = cell(nFiles, 1);
        for k = 1:nFiles
            
            %read stim file and get path to .wav file
            stimLinkFile = stimLinkFiles{k};
            fid = fopen(stimLinkFile);
            wavFile = strtrim(fgetl(fid));
            fclose(fid);
            
            stimFiles{k} = fullfile(allStimsDir, wavFile);            
        end
        srPairs.stimFiles = stimFiles;
        srPairs.respFiles = respFiles;
    end


