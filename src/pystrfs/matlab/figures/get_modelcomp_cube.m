function data = get_modelcomp_cube(cellName, preprocType, stimClass, fnames)

    if ~exist('fnames', 'var')
        fnames = {'directfit.%s.%s.%d.h5', 'directfit.nl_dists.%s.%s.%d.h5', ...
                  'directfit.nl_spline.%s.%s.%d.h5', 'threshgrad.%s.%s.0.50.%d.h5', ...
                  'threshgrad.%s.%s.exponential.0.50.%d.h5', ...
                  'threshgrad.%s.%s.logistic.0.50.%d.h5'};
    end

    parts = 0:9;
    
    rootDir = '/auto/k6/mschachter/pystrfs';
    
    switch preprocType
        case 'stft'
            preprocDesc = 'stft.nstd_6.fband_125';
        case 'lyons'
            [earQ, step] = get_best_lyons_params(cellName, stimClass);
            preprocDesc = sprintf('lyons.agc_1.earQ_%d.step_%0.2f', earQ, step);
        case 'surprise'
            preprocDesc = sprintf('surprise.dfw_3.dtw_3.dg_4.%s', stimClass);
    end
    
    strfs = cell(length(fnames), length(parts));
    strfMats = cell(length(fnames), 1);
    strfSizes = cell(length(fnames), 1);
    biases = zeros(length(fnames), length(parts));
    infos = zeros(length(fnames), length(parts), 2); %(:, : 1) = bound, (:, :, 2) = validation
    outputFiles = cell(length(fnames), length(parts));
    
    outputDir = fullfile(rootDir, 'units', cellName, 'output');
    
    h5 = h5utils();
    
    for k = 1:length(fnames)
       
        fname = fnames{k};
        strfSz = -1;
        for m = 1:length(parts)
            
            gnum = parts(m);
            outputFile = fullfile(outputDir, sprintf(fname, stimClass, preprocDesc, gnum));
                        
            if ~exist(outputFile, 'file')
                error('No such file: %s', outputFile);
            end
            
            fid = h5.open(outputFile);
            
            optMethod = h5.get_attr(fid, '/opt', 'method');
            if strcmp(optMethod, 'direct_fit') || strcmp(optMethod, 'threshgrad')
                strf = h5.get_ds(fid, '/model/weights');
                bias = h5.get_ds(fid, '/model/bias');            
            else
                strf = [];
                bias = 0.0;
            end
                        
            if strfSz == -1
                strfSz = size(strf);
            end
            
            infoBound = h5.get_ds(fid, '/model/performance/coherence/bound/info_mean');
            infoValid = h5.get_ds(fid, '/model/performance/coherence/validation/info_mean');
            
            strfs{k, m} = strf;
            biases(k, m) = bias;
            infos(k, m, 1) = infoBound;
            infos(k, m, 2) = infoValid;
            outputFiles{k, m} = outputFile;
            
            h5.close(fid);
        end
        
        strfSizes{k} = strfSz;
        %reshape strfs into a multidimensional array
        strfMat = zeros(strfSz(1), strfSz(2), length(parts));
        for m = 1:length(parts)        
            s = strfs{k, m};            
            strfMat(:, :, m) = s;
        end
        strfMats{k} = strfMat;
        
    end
    
    data = struct;
    data.fnames = fnames;
    data.strfs = strfs;
    data.strfMats = strfMats;
    data.biases = biases;
    data.infos = infos;
    data.strfSizes = strfSizes;
    data.outputFiles = outputFiles;
