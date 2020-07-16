function statdata = unit_to_statdata(unitFile, stimClass)

    unit = read_unit_file(unitFile);
    
    respData = unit.class_responses.(stimClass);
    
    statdata = struct;
    
    statdata.M = length(respData);
    statdata.N = 1;
    
    site = struct;
    site.label = {unitFile.description};
    site.time_scale = 1;
    site.time_resolution = 0.001;
    site.si_unit = 's';
    site.si_prefix = 1;
    
    statdata.sites(1) = site;
    
    
    
    