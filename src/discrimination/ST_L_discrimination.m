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


% Allocate file names for the temp results
done_mlnoise_name = cell(1,number_cells);
done_conspecific_name = cell(1,number_cells);

for nc=1:number_cells
    stimtype = 'flatrip';
    done_mlnoise_name{nc} = sprintf('%s_%s_%s_%s.mat' , fieldL_ST{nc,1}, brainregion, fieldL_ST{nc,2}, stimtype);
    stimtype = 'conspecific';
    done_conspecific_name{nc} = sprintf('%s_%s_%s_%s.mat' , fieldL_ST{nc,1}, brainregion, fieldL_ST{nc, 2}, stimtype);
end

not_done = 1;
n_loop = 0;
njobs = 0;
% Add jobs to queue until done
while not_done
    not_done = 0;
    for nc=1:number_cells
        stimtype = 'flatrip';
        if (exist(done_mlnoise_name{nc},'file') == 0)
            the_command = sprintf('cd(''/auto/fhome/fet/matlab/Neural Discrimination''); neural_discrimination_queue(''%s'', ''%s'', ''%s'', ''%s'')', fieldL_ST{nc, 1}, brainregion, fieldL_ST{nc, 2}, stimtype);
            the_comment = sprintf('ND%s%s%s%s',fieldL_ST{nc, 1}, brainregion, fieldL_ST{nc, 2}, stimtype);
            njobs = njobs +1;
            qid_submitted(njobs) = dbaddqueuemaster(the_command, the_comment);
            not_done = 1;
            if n_loop
                fprintf(1, 'Job %s not finished after loop %d\n', the_comment, n_loop);
            end
        end
        stimtype = 'conspecific';
        if (exist(done_conspecific_name{nc},'file') == 0)
            the_command = sprintf('cd(''/auto/fhome/fet/matlab/Neural Discrimination''); neural_discrimination_queue(''%s'', ''%s'', ''%s'', ''%s'')', fieldL_ST{nc, 1}, brainregion, fieldL_ST{nc, 2}, stimtype);
            the_comment = sprintf('ND%s%s%s%s',fieldL_ST{nc, 1}, brainregion, fieldL_ST{nc, 2}, stimtype);
            njobs = njobs +1;
            qid_submitted(njobs) = dbaddqueuemaster(the_command, the_comment);
            not_done = 1;
            if n_loop
                fprintf(1, 'Job %s not finished after loop %d\n', the_comment, n_loop);
            end
        end
    end
    pause(10*60); % pause for 10 minutes
    n_loop = n_loop +1;
end

% Clean the queue
for ij=1:njobs
    silent_dbdeletequeue(qid_submitted(ij));
end

% Print the results to a text file
% Allocate space for values in case we want to save them
nfiles_zfsong = zeros(1, number_cells);
stimrate_zfsong = zeros(1, number_cells);
stdstimrate_zfsong = zeros(1, number_cells);
backrate_zfsong = zeros(1, number_cells);
stdbackrate_zfsong = zeros(1, number_cells);
zscore_zfsong = zeros(1, number_cells);
stdz_zfsong = zeros(1, number_cells);
dprime_zfsong = zeros(1, number_cells);
stdd_zfsong = zeros(1, number_cells);
p_zfsong = zeros(1, number_cells);
rank_zfsong = zeros(1, number_cells);
info_zfsong =  zeros(1, number_cells);
noiseinfo_zfsong = zeros(1, number_cells);
rateinfo_zfsong = zeros(1, number_cells);

nfiles_mlnoise = zeros(1, number_cells);
stimrate_mlnoise = zeros(1, number_cells);
stdstimrate_mlnoise = zeros(1, number_cells);
backrate_mlnoise = zeros(1, number_cells);
stdbackrate_mlnoise = zeros(1, number_cells);
zscore_mlnoise = zeros(1, number_cells);
stdz_mlnoise = zeros(1, number_cells);
dprime_mlnoise = zeros(1, number_cells);
stdd_mlnoise = zeros(1, number_cells);
p_mlnoise = zeros(1, number_cells);
rank_mlnoise = zeros(1, number_cells);
info_mlnoise =  zeros(1, number_cells);
noiseinfo_mlnoise = zeros(1, number_cells);
rateinfo_mlnoise = zeros(1, number_cells);

info = 0; % info is used as a function as well - set it to zero to prevent strange Matlab behavior.
fid = fopen(fname, 'w');
fprintf(fid,'Bird\tCell\tNstim zf\tNstim ml\t');
fprintf(fid,'Stim Rate zf\tStim Rate ml\tSRstd zf\tSRstd ml\t');
fprintf(fid,'Back Rate zf\tBack Rate ml\tBRstd zf\tBRstd ml\t');
fprintf(fid,'Zscore zf\tZscore ml\tZstd zf\tZstd ml\t');
fprintf(fid,'Dprime zf\tDprime ml\tDstd zf\tDstd ml\t');
fprintf(fid,'P zf\tP ml\tRank zf\tRank ml\t');
fprintf(fid,'Info zf\tInfo ml\tInoise zf\tInoise ml\tIrate zf\tIrate ml\n');

for nc=1:number_cells
    
    load(done_conspecific_name{nc});
    nfiles_zfsong(nc) = nfiles;
    stimrate_zfsong(nc) = stimrate;
    stdstimrate_zfsong(nc) = stdstimrate;
    backrate_zfsong(nc) = backrate; 
    stdbackrate_zfsong(nc) = stdbackrate;
    zscore_zfsong(nc) = avgzscore;
    stdz_zfsong(nc) = stdzscore;
    dprime_zfsong(nc) = avgdprime;
    stdd_zfsong(nc) = stddprime;
    p_zfsong(nc) = percorrect;
    rank_zfsong(nc) = avgrank;
    info_zfsong(nc) = info;
    noiseinfo_zfsong(nc) = noiseinfo;
    rateinfo_zfsong(nc) = rateinfo;

    load(done_mlnoise_name{nc});
    nfiles_mlnoise(nc) = nfiles;
    stimrate_mlnoise(nc) = stimrate;
    stdstimrate_mlnoise(nc) = stdstimrate;
    backrate_mlnoise(nc) = backrate; 
    stdbackrate_mlnoise(nc) = stdbackrate;
    zscore_mlnoise(nc) = avgzscore;
    stdz_mlnoise(nc) = stdzscore;
    dprime_mlnoise(nc) = avgdprime;
    stdd_mlnoise(nc) = stddprime;
    p_mlnoise(nc) = percorrect;
    rank_mlnoise(nc) = avgrank;
    info_mlnoise(nc) = info;
    noiseinfo_mlnoise(nc) = noiseinfo;
    rateinfo_mlnoise(nc) = rateinfo;
    
    fprintf(fid,'%s\t%s\t%d\t%d\t', fieldL_ST{nc, 1}, fieldL_ST{nc, 2}, nfiles_zfsong(nc), nfiles_mlnoise(nc));    
    fprintf(fid,'%f\t%f\t%f\t%f\t', stimrate_zfsong(nc), stimrate_mlnoise(nc), stdstimrate_zfsong(nc), stdstimrate_mlnoise(nc) );
    fprintf(fid,'%f\t%f\t%f\t%f\t', backrate_zfsong(nc), backrate_mlnoise(nc), stdbackrate_zfsong(nc), stdbackrate_mlnoise(nc) );
    fprintf(fid,'%f\t%f\t%f\t%f\t', zscore_zfsong(nc), zscore_mlnoise(nc), stdz_zfsong(nc), stdz_mlnoise(nc) );
    fprintf(fid,'%f\t%f\t%f\t%f\t', dprime_zfsong(nc), dprime_mlnoise(nc), stdd_zfsong(nc), stdd_mlnoise(nc) );
    fprintf(fid,'%f\t%f\t%f\t%f\t', p_zfsong(nc), p_mlnoise(nc), rank_zfsong(nc), rank_mlnoise(nc) );
    fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n', info_zfsong(nc), info_mlnoise(nc), noiseinfo_zfsong(nc), noiseinfo_mlnoise(nc), ...
                                            rateinfo_zfsong(nc), rateinfo_mlnoise(nc));

end

fclose(fid);
