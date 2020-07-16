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
%   Author: Mike Schachter (mike.schachter@gmail.com)
%
function jobIds = slurm_sbatch(cmds, jobParams, matlabRunner)

    DEBUG = isfield(jobParams, 'debug');

    if nargin < 3
        matlabRunner = 'matlabbg';
    end
    
    if ~iscell(cmds)
        cmds = {cmds};
    end
       
    %% set defaults
    if ~isfield(jobParams, 'partition')
        jobParams.partition = 'all';
    end
    if ~isfield(jobParams, 'cpus')
        jobParams.cpus = 1;
    end
    
    %% construct list of commands to pass to sbatch
    sbatchCmds = {};
    if isfield(jobParams, 'partition')
        sbatchCmds{length(sbatchCmds)+1} = '-p';
        sbatchCmds{length(sbatchCmds)+1} = jobParams.partition;
    end
    if isfield(jobParams, 'depends')
        sbatchCmds{length(sbatchCmds)+1} = '-d';
	dstr = sprintf('afterok:%s', jobParams.depends);
        sbatchCmds{length(sbatchCmds)+1} = dstr;
    end
    if isfield(jobParams, 'out')
        sbatchCmds{length(sbatchCmds)+1} = '-o';
        sbatchCmds{length(sbatchCmds)+1} = jobParams.out;
    end
    if isfield(jobParams, 'err')
        sbatchCmds{length(sbatchCmds)+1} = '-e';
        sbatchCmds{length(sbatchCmds)+1} = jobParams.err;
    end
    if isfield(jobParams, 'nodes')
        sbatchCmds{length(sbatchCmds)+1} = '-N';
        sbatchCmds{length(sbatchCmds)+1} = num2str(jobParams.nodes);
    end
    if isfield(jobParams, 'qos')
        sbatchCmds{length(sbatchCmds)+1} = sprintf('--qos=%d', jobParams.qos);
    end
    if isfield(jobParams, 'cpus')
        sbatchCmds{length(sbatchCmds)+1} = '-c';
        sbatchCmds{length(sbatchCmds)+1} = num2str(jobParams.cpus);
    end
    if isfield(jobParams, 'nodelist')
        sbatchCmds{length(sbatchCmds)+1} = '-w';
        sbatchCmds{length(sbatchCmds)+1} = num2str(jobParams.nodelist);
    end
    
    
    sbatchRunCmd = 'sbatch';
    for k = 1:length(sbatchCmds)
        sbatchRunCmd = [sbatchRunCmd ' ' sbatchCmds{k}];
    end
    
    jobIds = ones(1, length(cmds)) * -1;
    
    for k = 1:length(cmds)
    
       %% write out executable command in script
       tempOut = tempname();
       fid = fopen(tempOut, 'w');
       fprintf(fid, '#!/bin/sh\n');
       fprintf(fid, '#SBATCH\n');
       
       mcmd = cmds{k};
       mcmd = strrep(mcmd, '\n', '');
       
       mThreadStr = sprintf('maxNumCompThreads(%d);', jobParams.cpus);
       
       mcmd = ['"' mcmd '"'];
       if DEBUG
	 fprintf('tempOut=%s\n', tempOut);
	 fprintf('mcmd=%s\n', mcmd);
       end
       fprintf(fid, '%s %s\n', matlabRunner, mcmd);
       fclose(fid);
       
       %% construct final sbatch command
       sbatchFinalCmd = [sbatchRunCmd ' ' tempOut];
       
       %% execute sbatch
       [status, result] = system(sbatchFinalCmd);
       if status ~= 0
           error('slurm_srun: problem executing command: %s', result);
       end 
       
       %% get job ID from result
       mat = regexp(result, 'job\s\d*', 'match');
       if ~isempty(mat)
           jobIds(k) = str2num(mat{1}(4:end));
       end
        
    end
       
