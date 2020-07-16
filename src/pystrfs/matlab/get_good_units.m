function units = get_good_units(units_and_regions_file)

    units_dir = '/auto/k6/mschachter/pystrfs/units';

    if nargin < 1 
        units_and_regions_file = fullfile(units_dir, 'good_units_and_regions.csv');
    end
    
    f = fopen(units_and_regions_file);    
    fdata = textscan(f, '%s %s', 'delimiter', ',');
    fclose(f);
    
    ncells = length(fdata{1});
    
    units = cell(ncells, 1);
    
    for k = 1:ncells
        
        unit_name = fdata{1}{k};
        unit_path = fullfile(units_dir, unit_name, 'unit.h5');
        region = fdata{2}{k};
        
        unit_obj = read_unit_file(unit_path);
        unit_obj.region = region;        
        
        units{k} = unit_obj;        
    end