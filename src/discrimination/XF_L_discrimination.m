% Runs the neural discrimination for cross-fostered L
brainregion = 'XFosteredZF_L';
fname = 'xf_fieldL.mat';

% All the data files pasted from the Excel file
fieldL_xf = { ...
{'ww1414', '2_B'}, ...
{'ww1414', '3_B'}, ...
{'ww1414', '4_B'}, ...
{'ww1414', '5_B'}, ...
{'ww1414', '6_B'}, ...
{'ww1414', '8_B'}, ...
{'ww1414', '9_B'}, ...
{'hpidg0712', '4_B'}, ...
{'hpidg0712', '5_B'}, ...
{'lblb0303', '3_B'}, ...
{'lblb0303', '4_B'}, ...
{'lblb0303', '5_B'}, ...
{'lblb0303', '6_B'}, ...
{'lblb0303', '7_B'}, ...
{'lblb0303', '8_B'}, ...
{'lblb0303', '9_B'}, ...
{'lblb0303', '10_B'}, ...
{'lblb0303', '11_B'}, ...
{'hpw0411', '4_B'}, ...
{'hpw0411', '6_B'}, ...
{'hpw0411', '7_B'}, ...
{'hpw0411', '8_B'}, ...
{'hpw0411', '9_B'}, ...
{'hpw0411', '10_B'}, ...
{'hpw0411', '11_B'}, ...
{'lblb0808', '3_B'}, ...
{'lblb0808', '4_B'}, ...
{'lblb0808', '5_B'}, ...
{'lblb0808', '6_B'}, ...
{'lblb0808', '7_B'}, ...
{'ww1212', '2_B'}, ...
{'ww1212', '3_B'}, ...
{'ww1212', '6_B'}, ...
{'ww1212', '12_B'}, ...
{'ww1212', '13_B'}, ...
};

% Allocate space for output values
number_cells = length(fieldL_xf);


% Allocate file names for the temp results
%done_bengalese_name = cell(1,number_cells);
done_conspecific_name = cell(1,number_cells);
%done_flatrip_name = cell(1,number_cells);

for nc=1:number_cells
    %stimtype = 'bengalese';
    %done_bengalese_name{nc} = sprintf('%s_%s_%s_%s.mat' , fieldL_xf{nc}{1}, brainregion, fieldL_xf{nc}{2}, stimtype);
    stimtype = 'conspecific';
    done_conspecific_name{nc} = sprintf('%s_%s_%s_%s.mat' , fieldL_xf{nc}{1}, brainregion, fieldL_xf{nc}{2}, stimtype);
    %stimtype = 'flatrip';
    %done_flatrip_name{nc} = sprintf('%s_%s_%s_%s.mat' , fieldL_xf{nc}{1}, brainregion, fieldL_xf{nc}{2}, stimtype);

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
        %    the_command = sprintf('cd(''/auto/fhome/smunro/matlab/Neural Discrimination''); neural_discrimination_queue(''%s'', ''%s'', ''%s'', ''%s'')', fieldL_xf{nc}{1}, brainregion, fieldL_xf{nc}{2}, stimtype);
        %    the_comment = sprintf('ND%s%s%s%s',fieldL_xf{nc}{1}, brainregion, fieldL_xf{nc}{2}, stimtype);
        %    njobs = njobs +1;
        %    qid_submitted(njobs) = dbaddqueuemaster(the_command, the_comment);
        %    not_done = 1;
        %    if n_loop
        %        fprintf(1, 'Job %s not finished after loop %d\n', the_comment, n_loop);
        %    end
        %end
        stimtype = 'conspecific';
        if (exist(done_conspecific_name{nc},'file') == 0)
            the_command = sprintf('cd(''/auto/fhome/smunro/matlab/Neural Discrimination''); neural_discrimination_queue(''%s'', ''%s'', ''%s'', ''%s'')', fieldL_xf{nc}{1}, brainregion, fieldL_xf{nc}{2}, stimtype);
            the_comment = sprintf('ND%s%s%s%s',fieldL_xf{nc}{1}, brainregion, fieldL_xf{nc}{2}, stimtype);
            njobs = njobs +1;
            qid_submitted(njobs) = dbaddqueuemaster(the_command, the_comment);
            not_done = 1;
            if n_loop
                fprintf(1, 'Job %s not finished after loop %d\n', the_comment, n_loop);
            end
        end
        %stimtype = 'flatrip';
        %if (exist(done_flatrip_name{nc},'file') == 0)
        %    the_command = sprintf('cd(''/auto/fhome/smunro/matlab/Neural Discrimination''); neural_discrimination_queue(''%s'', ''%s'', ''%s'', ''%s'')', fieldL_xf{nc}{1}, brainregion, fieldL_xf{nc}{2}, stimtype);
        %    the_comment = sprintf('ND%s%s%s%s',fieldL_xf{nc}{1}, brainregion, fieldL_xf{nc}{2}, stimtype);
        %    njobs = njobs +1;
        %    qid_submitted(njobs) = dbaddqueuemaster(the_command, the_comment);
        %    not_done = 1;
        %    if n_loop
        %        fprintf(1, 'Job %s not finished after loop %d\n', the_comment, n_loop);
        %    end
        %end

    end
    if ( not_done )
    pause(30*60); % pause for 30 minutes
    end
    n_loop = n_loop +1;
end

% Clean the queue
for ij=1:njobs
    silent_dbdeletequeue(qid_submitted(ij));
end

% Print the results to a text file
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
save(fname, 'nstim', 'stim_type', 'fieldL_xf', 'all_outputs');