function unit=read_unit_sectionsdata(h5Path)
    unit=read_unit_h5file(h5Path);
   
    h5 = h5utils();
    fid = h5.open(h5Path);
    
    crList = h5.get_subgroups(fid, '/');
    if ~ismember('extra_info', crList)
       return;
    end
    for k = 1:length(unit.classes)
        protocol=unit.classes{k};
        ciFatherList = h5.get_subgroups(fid, '/extra_info/');
        if ~ismember(protocol, ciFatherList)
            return
        end
        ciPath = sprintf('/extra_info/%s', protocol);
        ciList = h5.get_subgroups(fid, ciPath);
        if ~ismember('sections_data', ciList)
            return;
        end
        spath = sprintf('%s/sections_data', ciPath);

        cstruct.num_signif_sections_Bonferroni = h5.get_ds(fid, [spath '/num_signif_sections_Bonferroni']);
        cstruct.num_signif_sections_BonferroniDeg = h5.get_ds(fid, [spath '/num_signif_sections_BonferroniDeg']);
        cstruct.sections_pvalues = h5.get_ds(fid, [spath '/sections_pvalues']);
        unit.sections_data.(protocol) = cstruct;
    end
    