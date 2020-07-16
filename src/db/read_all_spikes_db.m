function [nfiles, spike_times, stimlen, stimname] = read_all_spikes_db(...
						classresp_id,respsamprate_Hz)
						
% Reads all the spike arrival times in a folder in the data bases as well as the lenght of
% the stimulus and the name of the stimulus.  The length of the stimulus is required to make sense of
% the data

% Default sampling rate is 1kHz
if nargin < 2
	respsamprate_Hz = 1000;
end

query_url = 'http://app.fet.berkeley.edu/query';
represp_ids = urlread(sprintf('%s/classresponse/%d/repeated.id',...
							  query_url,classresp_id));
represp_ids = sscanf(represp_ids,'%d');

if represp_ids(1) > 0

    % Read the list of spike data files in the folder and the stim length
    nfiles = length(represp_ids);

    % Read all the spike arrival times for each stim
    spike_times = cell(1, nfiles);
    stimname = cell(1, nfiles);
    stimlen = zeros(1, nfiles);

    for nf=1:nfiles

        % Get the database id number of the response
        idno = represp_ids(nf);

		% Get single trial presentation ids for this stim
		singleresp_ids = urlread(sprintf('%s/repeatedresponse/%d/singles.id',...
							  			 query_url,idno));
		singleresp_ids = sscanf(singleresp_ids,'%d');
		ntrials = length(singleresp_ids);

        % Read the spike arrival data
        for it=1:ntrials
            tline = urlread(sprintf('%s/singleresponse/%d/spiketimes_string',...
							  			 query_url,singleresp_ids(it)));
            if isempty(tline)
                spike_times{nf}{it} = NaN;
            else
                spike_times{nf}{it} = respsamprate_Hz * sscanf(tline,'%f');
            end
        end

        % Read the stimulus file to get its code and length
        stim_name = urlread(sprintf('%s/repeatedresponse/%d/presentation.stim.file.name',...
        							query_url,idno));
        							
        % Find only last part of stim_name
        remain = stim_name;
        while remain
            [token remain] = strtok(remain, '/');
        end
        stim_name = token;
        stimname{nf} = stim_name;
        
        stim_dur = urlread(sprintf('%s/repeatedresponse/%d/presentation.stim.file.duration',...
        							query_url,idno));
        stimlen(nf) = round(sscanf(stim_dur,'%f') * respsamprate_Hz);

    end
else
    fprintf(1, 'No stim presentations for class presentation %d',classpres_id);
    nfiles = 0;
end

if (nfiles == 0 )
    clear stimlen spike_times;
    stimlen = [];
    spike_times = [];
end