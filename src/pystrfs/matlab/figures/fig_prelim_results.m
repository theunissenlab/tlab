function fig_prelim_results(unit, preproc, model, thresh)
    
    preprocName = '';
    switch preproc
        case 'stft'
            preprocName = 'stft.nstd_6.fband_125';
        case 'lyons'
            preprocName = 'lyons.agc_1.earQ_8.step_0.50';
        case 'surprise'
            preprocName = 'surprise.dfw_3.dtw_3.dg_4.Con';
    end
    
    outputDir = sprintf('/auto/k6/mschachter/pystrfs/units/%s/output', unit);
    fname = sprintf('threshgrad.%s.Con.%s.thresh_%0.2f.0.h5', model, preprocName, thresh);
    outputFile = fullfile(outputDir, fname);
    
    m = get_model_data(outputFile);
        
    figure; hold on;
    imagesc(1:size(m.model.strf, 2), m.data.f, m.model.strf);
    cval = max(abs(m.model.strf(:)));
    caxis([-cval cval]);
    title(sprintf('STRF: bias=%0.4f', m.model.bias));
    axis tight;
    
    figure; hold on;
    xmin = min(m.model.response);
    xmax = max(m.model.response);
    xinc = (xmax - xmin) / 200;
    x = xmin:xinc:xmax;
    y = m.model.outputNL(x);
    plot(x, y, 'k-', 'LineWidth', 3); axis tight;
    axis([xmin xmax 0 max(m.data.wholeResp)]);
    title('Output NL');
    
    figure; hold on;
    plot(m.data.wholeResp, 'k-', 'LineWidth', 3);
    plot(m.model.response, 'r-', 'LineWidth', 2);
    axis([1 length(m.model.response) 0 max(m.data.wholeResp)]);
    legend('PSTH', 'Model');
    
    cb = m.coherence.bound;
    cm = m.coherence.validation;
    
    figure; hold on;
    plot(cb.f, cb.upper, 'r-', 'LineWidth', 3);
    plot(cb.f, cb.mean, 'k-', 'LineWidth', 3);
    plot(cb.f, cb.lower, 'b-', 'LineWidth', 3);
    
    plot(cm.f, cm.upper, 'r--', 'LineWidth', 3);
    plot(cm.f, cm.mean, 'k--', 'LineWidth', 3);
    plot(cm.f, cm.lower, 'b--', 'LineWidth', 3);    
    legend('Bound', 'Model');
    axis tight;
    