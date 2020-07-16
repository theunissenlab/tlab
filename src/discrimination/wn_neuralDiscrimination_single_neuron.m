% Runs the neural discrimination for wn-reared birds and neural recordings
% from field L
brainregion = 'wnL'; %NEED TO ADD THIS DATABASE TO dataCD_Feb09 directory in pgill's account
fname = 'wnreared_L_conMLnoise_info.mat';

%These are names of cells that I had "good" STRFs for wn birds

wn_allunits = { ...
    {'wnmale1', '5' }, ... %wnmale1 is 90 days old
    {'wnmale1', '6' }, ... %wnmale1 is 90 days old
    {'wnmale4', '9' }, ... 
    {'wnmale6', '8' }, ... 
    {'wnmale6', '9' }, ... 
    {'wnmale6', '10' }, ... 
    {'wnmale6', '11' }, ... 
    {'wnmale7', '11' }, ... 
    {'wnmale7', '15' }, ... 
    {'wnmale8', '1' }, ... 
    {'wnmale8', '2' }, ... 
    {'wnmale8', '3' }, ... 
    {'wnfemale1', '6' }, ... 
    {'wnfemale1', '7' }, ... 
    {'wnfemale1', '8' }, ... 
    {'wnfemale1', '9' }, ... 
    {'wnfemale4', '12' }, ... 
    {'wnfemale4', '13' }, ... 
    {'wnfemale5', '5' }, ... 
    {'wnfemale5', '7' }, ... 
    {'wnfemale5', '8' }, ... 
    {'wnfemale5', '9' }, ... 
    {'wnfemale5', '12' }, ... 
    {'wnfemale6', '5' }, ... 
    {'wnfemale6', '6' }, ... 
    {'wnfemale7', '13' }, ... 
    {'wnfemale8', '6' }};
    %{'wnfemale3', '5' }, ... %wnfemale3 is 8 months old
    %{'wnfemale3', '7' }, ... %wnfemale3 is 8 months old
    %{'wnfemale3', '8' }, ... %wnfemale3 is 8 months old
    %{'wnfemale3', '11' }}; %wnfemale3 is 8 months old

number_cells = length(wn_allunits);


% Allocate file names for the temp results
%done_bengalese_name = cell(1,number_cells);
done_con_name = cell(1,number_cells);
done_mlnoise_name = cell(1,number_cells);

for nc=1:number_cells
    %stimtype = 'bengalese';
    %done_bengalese_name{nc} = sprintf('%s_%s_%s_%s.mat' , fieldL_normal{nc}{1}, brainregion, fieldL_normal{nc}{2}, stimtype);
    stimtype = 'conspecific';
    done_con_name{nc} = sprintf('%s_%s_%s_%s.mat' , wn_allunits{nc}{1}, brainregion, wn_allunits{nc}{2}, stimtype);
    stimtype = 'flatrip';
    done_mlnoise_name{nc} = sprintf('%s_%s_%s_%s.mat' , wn_allunits{nc}{1}, brainregion, wn_allunits{nc}{2}, stimtype);
end

not_done = 1;
n_loop = 0;
njobs = 0;
% Add jobs to queue until done
while not_done
    not_done = 0;
    for nc=1:number_cells
        %stimtype = 'bengalese';
        %if (exist(done_bengalese_name{nc},'file') == 0)
        %    the_command = sprintf('cd(''/auto/fhome/fet/matlab/Neural Discrimination''); neural_discrimination_queue(''%s'', ''%s'', ''%s'', ''%s'')', fieldL_normal{nc}{1}, brainregion, fieldL_normal{nc}{2}, stimtype);
        %    the_comment = sprintf('ND%s%s%s%s',fieldL_normal{nc}{1}, brainregion, fieldL_normal{nc}{2}, stimtype);
        %    njobs = njobs +1;
        %    qid_submitted(njobs) = dbaddqueuemaster(the_command, the_comment);
        %    not_done = 1;
        %    if n_loop
        %        fprintf(1, 'Job %s not finished after loop %d\n', the_comment, n_loop);
        %    end
        %end
        stimtype = 'conspecific';
        if (exist(done_con_name{nc},'file') == 0)
            the_command = sprintf('cd(''/auto/fdata/noopur/wn/infoanalysis''); neural_discrimination_queue(''%s'', ''%s'', ''%s'', ''%s'')', wn_allunits{nc}{1}, brainregion, wn_allunits{nc}{2}, stimtype);
            the_comment = sprintf('ND%s%s%s%s',wn_allunits{nc}{1}, brainregion, wn_allunits{nc}{2}, stimtype);
            njobs = njobs +1;
            qid_submitted(njobs) = dbaddqueuemaster(the_command, the_comment);
            not_done = 1;
            if n_loop
                fprintf(1, 'Job %s not finished after loop %d\n', the_comment, n_loop);
            end
        end
        stimtype = 'flatrip';
        if (exist(done_mlnoise_name{nc},'file') == 0)
           the_command = sprintf('cd(''/auto/fdata/noopur/wn/infoanalysis''); neural_discrimination_queue(''%s'', ''%s'', ''%s'', ''%s'')', wn_allunits{nc}{1}, brainregion, wn_allunits{nc}{2}, stimtype);
           the_comment = sprintf('ND%s%s%s%s',wn_allunits{nc}{1}, brainregion, wn_allunits{nc}{2}, stimtype);
           njobs = njobs +1;
           qid_submitted(njobs) = dbaddqueuemaster(the_command, the_comment);
           not_done = 1;
           if n_loop
               fprintf(1, 'Job %s not finished after loop %d\n', the_comment, n_loop);
           end
        end
        end
    if not_done
        pause(30*60); % pause for 25 minutes
    end
    n_loop = n_loop +1;
end

% Clean the queue
%for ij=1:njobs
    %silent_dbdeletequeue(qid_submitted(ij));
%end

% Make a single data base
% Print the results to a text file
%nstim = 3;
nstim = 2;
%stim_type = {'bengalese', 'conspecific', 'flatrip'};
stim_type ={'conspecific', 'flatrip'};
all_outputs = cell(number_cells, nstim);

for nc=1:number_cells
    %load(done_bengalese_name{nc});
    %all_outputs{nc,1} = outputs;
    
    load(done_con_name{nc});
    all_outputs{nc,1} = outputs;
    
    load(done_mlnoise_name{nc});
    all_outputs{nc,2} = outputs;
end

cd /auto/fdata/noopur/wn/infoanalysis;
save(fname, 'nstim', 'stim_type', 'wn_allunits', 'all_outputs');