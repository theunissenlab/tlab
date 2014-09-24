function unames = get_good_unit_names(units_and_regions_file)

    units_dir = '/auto/k6/mschachter/pystrfs/units';

    if nargin < 1 
        units_and_regions_file = fullfile(units_dir, 'good_units_and_regions.csv');
    end
    
    f = fopen(units_and_regions_file);    
    fdata = textscan(f, '%s %s', 'delimiter', ',');
    fclose(f);
    
    ncells = length(fdata{1});
    
    unames = cell(1, ncells);
    for k = 1:ncells        
        unames{k} = fdata{1}{k};
    end
    