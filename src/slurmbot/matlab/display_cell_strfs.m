function display_cell_strfs(strfData)

    ptypes = {'stft', 'lyons', 'surprise'};
    stypes = {'conspecific', 'flatrip'};
    
    ncols = length(stypes);
    npreprocs = length(ptypes);

    figure('name', sprintf('STRFs: %s', strfData.cellName));
    for k = 1:length(ptypes)
        for m = 1:length(stypes)
            
            ptype = ptypes{k};
            stype = stypes{m};       
            indx = m + (k-1)*ncols;
            subplot(npreprocs, ncols, indx); hold on;
            
            strf = strfData.(ptype).(stype).strf;
            b = strfData.(ptype).(stype).bias;
            mval = max(abs(strf(:)));
            
            imagesc(strf); axis tight;
            caxis([-mval mval]); colorbar;
            title(sprintf('%s | %s | b=%0.4f', ...
                  ptypes{k}, stypes{m}, b));
        end
    end
