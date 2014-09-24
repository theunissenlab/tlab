% Runs the neural discrimination for normal L
brainregion = 'maybe_L';
fname = 'ST_fieldL.txt';
databasefilename = '/auto/fdata/anne/sarah_thane_infostuff/fieldL/con/databasefieldLcon.txt';

% Open data base text file to find all files
fdb = fopen(databasefilename);
if ( fdb == -1 )
    fprintf(1, 'Error: Could not open file %s\n', databasefilename);
end

number_cells = 0;
while true
    tline = fgetl(fdb);
    if ~ischar(tline), break; end
    number_cells = number_cells + 1;
end
fieldL_ST = cell(number_cells, 2);
frewind(fdb);

for it=1:number_cells
    tline = fgetl(fdb);
    if isempty(tline)
        fprintf(1, 'Error: Found blank line on database\n');
    else
       [fieldL_ST{it,1} remain] = strtok(tline);
       fieldL_ST{it,2} = strtok(remain);
    end
end
fclose(fdb);


% Check to see if all files exists...
% First conspecific

for nc=1:number_cells
    stimtype = 'conspecific';
    [nfiles spike_times stim_len] = read_all_spikes(fieldL_ST{nc, 1}, brainregion, fieldL_ST{nc, 2}, stimtype);
% Then mlnoise
    stimtype = 'flatrip';
    [nfiles spike_times stim_len] = read_all_spikes(fieldL_ST{nc, 1}, brainregion, fieldL_ST{nc, 2}, stimtype);
end
