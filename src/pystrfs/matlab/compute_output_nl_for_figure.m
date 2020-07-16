function compute_output_nl_for_figure(responseFile)

    fprintf('Computing output NL data for %s\n', responseFile);
    m = get_model_data(responseFile);
    
    xmin = min(m.model.response);
    xmax = max(m.model.response);
    xinc = (xmax - xmin) / 200;
    x = xmin:xinc:xmax;
    y = m.model.outputNL(x);
    if isempty(x)
        x = 0;
    end
    if isempty(y)
        y = 0;
    end

    h5 = h5utils();
    fid = h5.open(responseFile);    
    h5.set_ds(fid, '/model/output_nl', 'domain', x);
    h5.set_ds(fid, '/model/output_nl', 'range', y);
    h5.close(fid);