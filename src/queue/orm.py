from datetime import datetime
from sqlalchemy import Table, Column
from sqlalchemy.orm import mapper
from .meta import metadata, Session

class Job(object):
    
    def __init__(self,parmstring='',progname='echo 1',
                 rundataid=0,priority=1,allowqueuemaster=1,user='DBUSER',
                 note='',waitid=0,memrequest=0):
        
        self.parmstring = parmstring
        self.rundataid = rundataid
        self.progname = progname
        self.priority = priority
        self.queuedate = datetime.now()
        self.allowqueuemaster = allowqueuemaster
        self.user = user
        self.note = note
        self.waitid = waitid
        self.memoryrequest = memrequest
    
tqueue = Table('tQueue',metadata,autoload=True,schema='cell')

mapper(Job,tqueue)

jobs = Session.query(Job)