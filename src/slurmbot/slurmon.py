import time
import logging

from slurmtools import *


#LOG_FILENAME = '/auto/k1/slurm/slurmon.log'
LOG_FILENAME = '/auto/fhome/mschachter/slurmon.log'
logging.basicConfig(filename=LOG_FILENAME, level=logging.DEBUG)


TOTAL_JOB_THRESH = 5
TOTAL_TIME_THRESH = 12*60*60

class Slurmon:
      
   pollInterval = 60.0

   def __init__(self):
      pass


   def run(self):
      
      while True:

         userStats = get_user_level_stats()
         jobsToDowngrade = []
         
         for (uname,ustats) in userStats.iteritems():
            
            jobInfos = slurm_squeue(username=uname)                  
            for ji in jobInfos:
               tstr = ji['TIME']
               tsecs = parse_time(tstr)
               ncpuTotal = ustats['TOTAL_CPUS']
               if tsecs > TOTAL_TIME_THRESH and ncpuTotal > TOTAL_JOB_THRESH:
                  jobsToDowngrade.append(ji['ID'])

         print jobsToDowngrade

         time.sleep(self.pollInterval)
         



if __name__ == '__main__':
   
   slurmon = Slurmon()
   slurmon.run()




