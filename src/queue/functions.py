import os
from . import meta
from .orm import Job, jobs

def command(*args,**kwargs):
    command = ' '.join(str(a) for a in args)
    dbaddqueuemaster(command,**kwargs)

def script(scriptname,*args,**kwargs):
    command('python',scriptname,*args,**kwargs)

def dbaddqueuemaster(command,**kwargs):
    '''
    function queueidx = dbaddqueuemaster(parmstring,note='',priority=1,
                                         rundataid=0,memrequest=0);

    add a job to the queue
    
    parameters: 
    parmstring - string containing command(s) for python
               
    note - note to appear for this job on the queue monitor
    priority - (0-5) higher values cause jobs to be executed sooner than
             jobs with lower values (among your OWN jobs)
    
    options -    other options (currently memory request only)
               'memory': request amount of memory / proc (in GB).
                       The job is assigned only to the machines
                       that have physical memory(/proc) more than
                       requested.
    
    ex.
     dbaddqueuemaster('c=1;calcx;', 'calc1', 1, 'memory', 1.5);
     
    Ported from /auto/k1/queued/dbaddqueuemaster.m 9/13/10 RCM
    '''
    kwargs['parmstring'] = ''
    kwargs['progname'] = command
    kwargs.setdefault('note','')
    kwargs.setdefault('priority',1)
    kwargs.setdefault('rundataid',0)
    kwargs.setdefault('memrequest',0)
    kwargs.setdefault('waitid',0)
    kwargs.setdefault('user',os.environ['USER'])
    
    j = Job(**kwargs)
    meta.Session.add(j)
    meta.Session.commit()
    meta.Session.refresh(j)

    return j.id

def get_jobs(**criteria):
    '''
    query = get_jobs(**criteria)
    
    Returns an SQLAlchemy query object with the parameters specified.
    criteria keywords may be any of the fields on table tQueue, most likely
    one of [user,complete,priority]
    
    e.g. get_jobs(user='channing',complete=-1) will return a query of all jobs
    with user channing that are currently running (complete = -1).
    
    The SQLAlchemy query object is iterable and supports indexing, though
    for some reason minus indexing don't work. E.g., q[10] works, but q[-1]
    does not.
    
    The query object also has other useful methods, most importantly:
    
        one(): return exactly one object; raise an error if the query returns
        more or fewer than one. Usually used as a sanity check when you should
        only get one record, e.g. filter_by(id=938).one()
        
        count(): return an integer corresponding to the number of records. 
    '''
    return jobs.filter_by(**criteria)