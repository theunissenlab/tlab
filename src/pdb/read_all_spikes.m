function [nfiles, spike_times, stimlen, stimname] = read_all_spikes(birdname, brainregion, cellname, stimtype)
% Reads all the spike arrival times in a folder in the data bases as well as the lenght of
% the stimulus and the name of the stimulus.  The length of the stimulus is required to make sense of
% the data

savecd = pwd;

% The rootcd will be different on the Linux Cluster versus individual
% machines
rootcd = '/auto/fdata/pgill/dataCD_April09'; % Linux Cluster
%rootcd = 'C:\Users\Frederic\Documents\Data\dataCD_withXF'; % Frederic's Paris Laptop

datacd = strcat(rootcd, '/', brainregion, '/', birdname, '_', cellname, '/', stimtype);

if exist(datacd, 'file' )
    cd(datacd);

    % Read the list of spike data files in the folder and the stim length
    datname = dir('spike*');
    nfiles = length(datname);

    % Read all the spike arrival times for each stim
    spike_times = cell(1, nfiles);
    stimname = cell(1, nfiles);
    stimlen = zeros(1, nfiles);

    for nf=1:nfiles

        % Get the decimal id number of the file
        idno = sscanf(datname(nf).name,'spike%d');

        % Read the spike arrival data
        fid = fopen(datname(nf).name);
        if ( fid == -1 )
            fprintf(1, 'Error: could not open file name %s/%s\n', datacd, datname(nf).name);
            nfiles = 0;
            break;
        end
        nt = 0;
        while true
            tline = fgetl(fid);
            if ~ischar(tline), break; end
            nt = nt + 1;
        end

        spike_times{nf} = cell(1,nt);
        frewind(fid);

        for it=1:nt
            tline = fgetl(fid);
            if isempty(tline)
                spike_times{nf}{it} = NaN;
            else
                spike_times{nf}{it} = sscanf(tline, '%f');
            end
        end
        fclose(fid);

        % Read the stimulus file to get its code and length
        fid = fopen(sprintf('stim%d', idno));
        if ( fid == -1 )
            fprintf(1, 'Error: could not open file name %s/stim%d\n', datacd, idno);
            nfiles = 0;
            break;
        end
        stim_name = fgetl(fid);
        
        fclose(fid);

        % Find only last part of stim_name
        remain = stim_name;
        while remain
            [token remain] = strtok(remain, '/');
        end
        stim_name = token;
        stimname{nf} = stim_name;
        fid = fopen(sprintf('%s/all_stim_lengths/%s', rootcd, stim_name));
        if (fid == -1 )
            fprintf(1, 'Error: could not open file name %s/all_stim_lengths/%s\n', rootcd, stim_name);
            nfiles = 0;
            break;
        end
        stimlen(nf) = fscanf(fid,'%d');
        fclose(fid);

    end
else
    fprintf(1, 'Error: Missing Data directory %s\n', datacd);
    nfiles = 0;
end

if (nfiles == 0 )
    clear stimlen spike_times;
    stimlen = [];
    spike_times = [];
end

cd(savecd);
return
