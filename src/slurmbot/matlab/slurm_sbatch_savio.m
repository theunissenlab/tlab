%% Run a string of MATLAB commands using slurm
%
%   Input:
%       cmds: Either a single string of commands, or a
%           cell array of strings. If a cell array, each
%           element is run using it's own job.
%
%       jobParams: Parameters for the jobs:
%           .partition: the name of the partition to run on
%           .depends: a comma-separated list of job id dependencies
%           .out: file path to write stdout to
%           .err: file path to write stderr to
%           .nodes: the # of nodes this job will require
%           .qos: the quality-of-service for this job
%           .cpus: the # of CPUs required by this job
%           .nodelist: a comma-separated list of specific nodes to run on
%
%       matlabRunner: (optional) the script used to run matlab commands,
%           defaults to 'matlabbg'
%
%   Returns:
%       jobIds: a vector of job IDs, the same size as cmds
%
%   Author:  Julie Elie
%
function jobIds = slurm_sbatch_savio(cmds, jobParams, matlabRunner)

    DEBUG = isfield(jobParams, 'debug');

    if nargin < 3
        matlabRunner = '/global/home/users/jelie/CODE/bin/matlabbg';
        %matlabRunner = 'matlab';
    end
    
    if ~iscell(cmds)
        cmds = {cmds};
    end
       
    %% set defaults
    if ~isfield(jobParams, 'Name')
        jobParams.Name = 'TheunissenLabJob';
    end
    if ~isfield(jobParams, 'Partition')
        jobParams.Partition = 'savio';
    end
    if ~isfield(jobParams, 'Account')
        jobParams.Account = 'fc_birdpow';
    end
    if ~isfield(jobParams, 'qos')
        jobParams.Qos = 'savio_normal';
    end
    if ~isfield(jobParams, 'NTasks')
        jobParams.NTasks = 1;
    end
    if ~isfield(jobParams, 'CPU')
        jobParams.CPU = 20; % You want this number to be the number of...
                            %iterations in your parallel loops up to 20 or...
                            %24 per job (there are 20 cores per node=computer...
                            %on savio and 24 cores per node on savio2)
    end
    if ~isfield(jobParams, 'TimeLimit')
        jobParams.TimeLimit = '72:00:00';
    end
    
    
    
    jobIds = ones(1, length(cmds)) * -1;
    
    %% write out executable command in script
    % Name for the file
    if ~isfield(jobParams, 'JobName')
        [Pathstr,Name,Ext]=fileparts(cmds{1});
    else
        Name = jobParams.JobName;
    end
       if isfield(jobParams, 'Type')
           tempOut = sprintf('ExJob_%s_%s.txt', Name, jobParams.Type);
       else
           tempOut = ['ExJob' Name '.txt'];
       end
       fid = fopen(tempOut, 'w');
       fprintf(fid, '#!/bin/bash\n');
       fprintf(fid, '# Job name:\n#SBATCH --job-name=%s\n#\n',jobParams.Name);
       fprintf(fid, '# Partition:\n#SBATCH --partition=%s\n#\n',jobParams.Partition);
       fprintf(fid, '# Account:\n#SBATCH --account=%s\n#\n',jobParams.Account);
       fprintf(fid, '# QoS:\n#SBATCH --qos=%s\n#\n',jobParams.Qos);
       fprintf(fid, '# Processors:\n#SBATCH --ntasks=%d\n#\n',jobParams.NTasks);
       fprintf(fid, '# Core per processor:\n#SBATCH --cpus-per-task=%d\n#\n',jobParams.CPU);
       fprintf(fid, '# Wall clock limit:\n#SBATCH --time=%s\n#\n',jobParams.TimeLimit);
       fprintf(fid, '# Error file:\n#SBATCH --error=%s\n#\n',jobParams.err);
       fprintf(fid, '# Output file:\n#SBATCH --output=%s\n#\n',jobParams.out);
       
       mcmd = strrep(cmds{1}, '\n', '');

       mcmd = ['"' mcmd '"'];
       if DEBUG
	 fprintf('tempOut=%s\n', tempOut);
	 fprintf('mcmd=%s\n', mcmd);
       end
       %fprintf(fid, '## Run command\nmodule load matlab\n%s -nosplash -nodesktop -logfile %s < %s\necho "end of Batch script"\nexit\n',matlabRunner, jobParams.out, mcmd);
       fprintf(fid, '## Run command\nmodule load matlab\n%s %s\n',matlabRunner, mcmd);
       %fprintf(fid, '## Run command\n%s %s\n',matlabRunner, mcmd);
       fclose(fid);
       
%        %% construct final sbatch command
%        sbatchFinalCmd = tempOut;
%        
%        %% execute sbatch
%        [status, result] = system(sbatchFinalCmd);
%        if status ~= 0
%            error('slurm_srun: problem executing command: %s', result);
%        end 
       
%        %% get job ID from result
%        mat = regexp(result, 'job\s\d*', 'match');
%        if ~isempty(mat)
%            jobIds(k) = str2num(mat{1}(4:end));
%        end
        
    end
       
