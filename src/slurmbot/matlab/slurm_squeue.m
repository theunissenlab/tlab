function jobInfo = slurm_squeue()
    
    cmd = 'squeue -h -o "%N %P %Q %u %M %T %i"';
    
    [stat, result] = system(cmd);
    
    result
    lns = string_split(result, '\n')
    
    jobInfo = struct;

    
    
    
end
    
function sarr = string_split(str, delim)

    indxs = strfind(str, delim);
    if ~isempty(indxs)
        sarr = cell(length(indxs)+1, 1);

        startIndx = 1;
        for k = 1:length(indxs)        
           i = indxs(k);
           sarr{k} = str(startIndx:(i-1));
           startIndx = i + length(delim);        
        end
        sarr{end} = str( (indxs(end)+length(delim)):end );
    else
        sarr = {str};
    end
    
end