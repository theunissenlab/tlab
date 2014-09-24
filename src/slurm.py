from os import system

def batch(scriptname,*args,**kwargs):
    '''
    Run script as a slurm batch (i.e. queued)
    *args are the bind parameters for the script,
    and must all be strings.
    '''
    
    outfile = kwargs.pop('outfile','~/slurm/out%j.txt')
    errfile = kwargs.pop('errfile','~/slurm/err%j.txt')
    partition = kwargs.pop('partition','all')
    share = kwargs.pop('share',False)
    nice = kwargs.pop('nice',8)
    nodes = kwargs.pop('nodes',2)
    
    switches = "-o %s -e %s -p %s -c %d" % (outfile,errfile,partition,nodes)
    if share:
        switches += ' -s'

    if nice is not None:
        switches += ' --nice=%d' % nice
    
    for switch,value in kwargs.items():
        switches += ' --%s' % switch
        if value is not None:
            switches += '=%s' % value
    
    system("sbatch %s %s %s" % (switches,
                                scriptname,
                                ' '.join(str(a) for a in args)))