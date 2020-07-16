function create_best_models_preproc(bestFile, outputFile)

    fid = fopen(bestFile, 'r');
    data = textscan(fid,'%s%s%s%s%f%f', 'Delimiter', ',');
    fclose(fid);
    
    unitsDir = '/auto/k6/mschachter/pystrfs/units';
    
    pname.surprise = 'surprise.dfw_3.dtw_3.dg_4.Con';
    pname.lyons = 'lyons.agc_1.earQ_8.step_0.50';
    pname.rawspectrogram = 'stft.nstd_6.fband_125';
    
    ftemp = 'threshgrad.v5.%s.Con.%s.thresh_%0.2f.%%d.h5';
    
    nunits = length(data{1});
    
    for k = 1:nunits
       
        unit = data{1}{k};
        region = data{2}{k};
        preproc = data{3}{k};
        model = data{4}{k};
        thresh = data{5}(k);
        score = data{6}(k);
        
        preprocName = pname.(preproc);
        
        if strcmp(model, 'sepnl_spline')
            model = 'linear.nl_spline';
        elseif strcmp(model, 'sepnl_dists')
            model = 'linear.nl_dists';
        end
        
        fname = sprintf(ftemp, model, preprocName, thresh);
        fpath = fullfile(unitsDir, unit, 'output', fname);
        
        preproc_model_response(fpath, unit, 'best', outputFile);
    end