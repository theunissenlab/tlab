function stimMd5s = get_stim_md5s(unitFile, stimClass)

    if nargin < 2
        stimClass = 'Con';
    end

    unit = read_unit_file(unitFile);
    
    resps = unit.class_responses.(stimClass);
    nresps = length(resps);
    stimMd5s = cell(nresps, 1);
    for k = 1:nresps        
        stimMd5s{k} = resps{k}.stimMd5;        
    end
    