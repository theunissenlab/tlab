function [stimFiles, md5s] = read_stimdata_from_file(fileName)

    fid = fopen(fileName);
    
    tline = fgetl(fid);
    stimFiles = {};
    md5s = {};
    while ischar(tline)
              
        cindx = strfind(tline, ',');
        
        md5 = tline(1:(cindx-1));
        fname = tline((cindx+1):end);
        
        stimFiles{end+1} = fname;
        md5s{end+1} = md5;
        
        tline = fgetl(fid);
    end
