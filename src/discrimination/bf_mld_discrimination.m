% Runs the neural discrimination for bf mld
brainregion = 'bengalese_mld';
fname = 'bf_mld.mat';

% All the data files pasted from the Excel file
mld_bf = { ...
{'pi05', '2_A'}, ...
{'pi05', '3_A'}, ...
{'pi05', '4_A'}, ...
{'pi05', '5_A'}, ...
{'pi05', '6_A'}, ...
{'pi05', '7_A'}, ...
{'pi05', '8_A'}, ...
{'pi05', '9_A'}, ...
{'pi05', '10_A'}, ...
{'pi05', '11_A'}, ...
{'g08', '2_A'}, ...
{'g08', '3_A'}, ...
{'g08', '6_A'}, ...
{'g08', '7_A'}, ...
{'gg1717', '2_A'}, ...
{'gg1717', '3_A'}, ...
{'gg1717', '4_A'}, ...
{'gg1717', '5_A'}, ...
{'gg1717', '6_A'}, ...
{'gg1717', '7_A'}, ...
{'gg1717', '8_A'}, ...
{'gg1717', '9_A'}, ...
{'gg1717', '10_A'}, ...
{'gg1717', '11_A'}, ...
{'pipi0607', '3_A'}, ...
{'pipi0607', '4_A'}, ...
{'pipi0607', '5_A'}, ...
{'pipi0607', '6_A'}, ...
{'pipi0607', '7_A'}, ...
{'pipi0607', '8_A'}, ...
{'pipi0607', '9_A'}, ...
{'g06', '2_A'}, ...
{'g06', '3_A'}, ...
{'g06', '4_A'}, ...
{'g06', '5_A'}, ...
{'g06', '6_A'}, ...
{'g03', '2_A'}, ...
{'g03', '5_A'}, ...
{'g03', '8_A'}, ...
{'g03', '9_A'}, ...
{'g03', '12_A'}, ...
{'lblb0505', '2_A'}, ...
{'lblb0505', '3_A'}, ...
{'lblb0505', '4_A'}, ...
{'lblb0505', '5_A'}, ...
{'lblb0505', '7_A'}, ...
{'lblb0505', '8_A'}, ...
{'lblb0505', '9_A'}, ...
{'lblb0505', '10_A'}, ...
{'pi07', '2_A'}, ...
{'pi07', '3_A'}, ...
{'pi07', '4_A'}, ...
{'pi07', '5_A'}, ...
{'pi07', '6_A'}, ...
{'pi07', '7_A'}, ...
{'pi07', '8_A'}, ...
{'pi07', '9_A'}};

% Allocate space for output values
number_cells = length(mld_bf);



% Allocate file names for the temp results
%done_bengalese_name = cell(1,number_cells);
done_conspecific_name = cell(1,number_cells);
%done_flatrip_name = cell(1,number_cells);

for nc=1:number_cells
    %stimtype = 'bengalese';
    %done_bengalese_name{nc} = sprintf('%s_%s_%s_%s.mat' , mld_bf{nc}{1}, brainregion, mld_bf{nc}{2}, stimtype);
    stimtype = 'conspecific';
    done_conspecific_name{nc} = sprintf('%s_%s_%s_%s.mat' , mld_bf{nc}{1}, brainregion, mld_bf{nc}{2}, stimtype);
    %stimtype = 'flatrip';
    %done_flatrip_name{nc} = sprintf('%s_%s_%s_%s.mat' , mld_bf{nc}{1}, brainregion, mld_bf{nc}{2}, stimtype);
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
        %    the_command = sprintf('cd(''/auto/fhome/smunro/matlab/Neural Discrimination''); neural_discrimination_queue(''%s'', ''%s'', ''%s'', ''%s'')', mld_bf{nc}{1}, brainregion, mld_bf{nc}{2}, stimtype);
        %    the_comment = sprintf('ND%s%s%s%s',mld_bf{nc}{1}, brainregion, mld_bf{nc}{2}, stimtype);
        %    njobs = njobs +1;
        %    qid_submitted(njobs) = dbaddqueuemaster(the_command, the_comment);
        %    not_done = 1;
        %    if n_loop
        %        fprintf(1, 'Job %s not finished after loop %d\n', the_comment, n_loop);
        %    end
        %end
        stimtype = 'conspecific';
        if (exist(done_conspecific_name{nc},'file') == 0)
            the_command = sprintf('cd(''/auto/fhome/smunro/matlab/Neural Discrimination''); neural_discrimination_queue(''%s'', ''%s'', ''%s'', ''%s'')', mld_bf{nc}{1}, brainregion, mld_bf{nc}{2}, stimtype);
            the_comment = sprintf('ND%s%s%s%s',mld_bf{nc}{1}, brainregion, mld_bf{nc}{2}, stimtype);
            njobs = njobs +1;
            qid_submitted(njobs) = dbaddqueuemaster(the_command, the_comment);
            not_done = 1;
            if n_loop
                fprintf(1, 'Job %s not finished after loop %d\n', the_comment, n_loop);
            end
        end
        %stimtype = 'flatrip';
        %if (exist(done_flatrip_name{nc},'file') == 0)
        %    the_command = sprintf('cd(''/auto/fhome/smunro/matlab/Neural Discrimination''); neural_discrimination_queue(''%s'', ''%s'', ''%s'', ''%s'')', mld_bf{nc}{1}, brainregion, mld_bf{nc}{2}, stimtype);
        %    the_comment = sprintf('ND%s%s%s%s',mld_bf{nc}{1}, brainregion, mld_bf{nc}{2}, stimtype);
        %    njobs = njobs +1;
        %    qid_submitted(njobs) = dbaddqueuemaster(the_command, the_comment);
        %    not_done = 1;
        %    if n_loop
        %        fprintf(1, 'Job %s not finished after loop %d\n', the_comment, n_loop);
        %    end
        %end
    end
    if (not_done)
        pause(30*60); % pause for 30 minutes
    end
    n_loop = n_loop +1;
end

% Clean the queue
for ij=1:njobs
    silent_dbdeletequeue(qid_submitted(ij));
end


% Print the results to a mat file
%nstim = 3;
nstim = 1;
%stim_type = {'bengalese', 'conspecific', 'flatrip'};
stim_type = {'conspecific'};
all_outputs = cell(number_cells, nstim);

for nc=1:number_cells
    %load(done_bengalese_name{nc});
    %all_outputs{nc,1} = outputs;
    
    load(done_conspecific_name{nc});
    all_outputs{nc,2} = outputs;
    
    %load(done_flatrip_name{nc});
    %all_outputs{nc,3} = outputs;
end
save(fname, 'nstim', 'stim_type', 'mld_bf', 'all_outputs');

