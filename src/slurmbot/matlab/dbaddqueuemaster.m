function jobIds = dbaddqueuemaster(paramstring, note, priority, varargin)

    if nargin < 2
        note = '';
    end
    if nargin < 3
        priority = 1;
    end

    memRequested = 1;
    narg = 1;
    while narg <= length(varargin)
        switch varargin{narg}
            case 'memory'
              memRequested = varargin{narg+1};
              narg = narg + 1;
        end
    end
    nCpus = round(ceil(memRequested));
              
    jobParams = struct;
    jobParams.comment = note;
    jobParams.cpus = nCpus;
    
    if priority <= 0
        jobParams.partition = 'lowprio';        
    elseif priority == 1
        jobParams.partition = 'all';
    else
        jobParams.partition = 'highprio';
    end
        
    jobIds = slurm_sbatch(paramstring, jobParams);
    