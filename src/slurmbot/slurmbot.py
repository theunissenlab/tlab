import sys
import string
import os
import slurmtools
import time

class SlurmJob:

   def __init__(self):
      self.cmds = None
      self.params = None
      self.completionFile = None
      self.running = False
      self.startTimes = 0
      self.completed = False
      self.jobId = -1
      

class SlurmBot:

   def __init__(self):
      self.jobs = []
      self.pollInterval = 37.0  
      self.maxJobs = 10


   def add(self, cmdStrs, sbatchParams, completionFileName = None):
      sj = SlurmJob()
      sj.cmds = cmdStrs
      sj.params = sbatchParams
      sj.completionFile = completionFileName
      self.jobs.append(sj)


   def clear(self):
      del self.jobs[0:len(self.jobs)]
      

   def get_queued_jobs(self):      
      return filter(lambda sj: not sj.completed and not(sj.running), self.jobs)


   def get_running_jobs(self):      
      return filter(lambda sj: sj.running and not(sj.completed), self.jobs)


   def get_uncompleted_jobs(self):
      ujList = []
      for sj in self.jobs:
         if not sj.completed or not os.path.exists(sj.completionFile):
            ujList.append(sj)
      return ujList

   def mark_completed_jobs(self):
      for sj in self.jobs:
         if sj.completionFile != None and os.path.exists(sj.completionFile):
            sj.completed = True

   def update_running_jobs(self):

      jobInfo = slurmtools.slurm_squeue()
      for rj in self.get_running_jobs():
         jinfo = filter(lambda j: j['ID'] == rj.jobId, jobInfo)
         if len(jinfo) == 0:
            rj.running = False
            if rj.completionFile is None:
               print 'Marking job %d as completed (no completion file specified)' % rj.jobId
               rj.completed = True
            elif rj.completionFile != None and os.path.exists(rj.completionFile):
               print 'Marking job %d as completed' % rj.jobId
               rj.completed = True

        
   def run_and_wait(self, ignoreCompleted = False):

      if not ignoreCompleted:
         self.mark_completed_jobs()
         uj = self.get_uncompleted_jobs()
         nComp = len(self.jobs) - len(uj)
         print 'Marked %d jobs as already completed' % nComp
      
      self.update_running_jobs()
      queuedJobs = self.get_queued_jobs()
      runningJobs = self.get_running_jobs()

      while len(queuedJobs) > 0 or len(runningJobs) > 0:
         
         self.update_running_jobs()
         
         queuedJobs = self.get_queued_jobs()
         runningJobs = self.get_running_jobs()

         print '# of queued jobs: %d' % len(queuedJobs)
         print '# of running jobs: %d' % len(runningJobs)

         nJobsAvailable = min(len(queuedJobs), self.maxJobs - len(runningJobs))
         print '# of available job space to run: %d' % nJobsAvailable
         if nJobsAvailable > 0 and len(queuedJobs) > 0:

            for k in range(nJobsAvailable):
               #Run jobs
               sj = queuedJobs[k]
               jobId = slurmtools.slurm_sbatch(sj.cmds, **sj.params)
               sj.running = True
               sj.jobId = jobId
               if sj.startTimes == 0:
                  print 'Started job with id %d' % jobId
               else: 
                  print 'Restarted job with id %d' % jobId
               sj.startTimes += 1
         
         time.sleep(self.pollInterval)


   def run(self):
      
      for sj in self.get_queued_jobs():
         cmdStrs = sj.cmds
         slurmParams = sj.params
         jobId = slurmtools.slurm_sbatch(cmdStrs, **slurmParams)
         sj.jobId = jobId
         print 'Started job with id %d' % jobId
       

if __name__ == '__main__':

   sb = SlurmBot()
   sb.maxJobs = 5
   
   cmdStr = ['sleep', '10s']
   sparams = {'out': '/auto/fhome/mschachter/slurm_%j_out.txt',
              'err': '/auto/fhome/mschachter/slurm_%j_err.txt',
              'partition':'all'}

   print 'adding processes'
   for k in range(10):
      sb.add(cmdStr, sparams)
   
   print 'running and waiting'   
   sb.run_and_wait()

   print 'done'

