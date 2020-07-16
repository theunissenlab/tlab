
import sys
import time

from modelcomps import *
from run_threshgrad import *
from run_lars import *
from run_directfit import *


if __name__ == '__main__':

   
   if len(sys.argv) < 3:
      print 'Usage:'
      print '  python run_lars.py preproc <cell1,cell2,...>'      
      print '  python run_lars.py run <cell1,cell2,...>'
      print '  python run_lars.py missing <cell1,cell2,...>'
      exit()

   cmd = sys.argv[1]

   if cmd == 'preproc':      
      cellNamesStr = sys.argv[2]
      cellNames = cellNamesStr.split(',')
      
      slurmBot = slurmbot.SlurmBot()
      slurmBot.maxJobs = 10

      for cellName in cellNames:
         print 'Cleaning and preprocessing cell %s...' % cellName
         add_preprocs(cellName, slurmBot, True, 'strflab*')         
  
      slurmBot.run()

      exit()

   if cmd == 'run':
      
      cellNamesStr = sys.argv[2]
      cellNames = cellNamesStr.split(',')
      
      slurmBot = slurmbot.SlurmBot()
      slurmBot.maxJobs = 7;

      for cellName in cellNames:
         add_directfit_optimizations(cellName, 'linear', slurmBot)
         #add_directfit_optimizations(cellName, 'exponential', slurmBot)
         add_lars_optimizations(cellName, 'linear', slurmBot)
         #add_lars_optimizations(cellName, 'exponential', slurmBot)
         add_threshgrad_optimizations(cellName, 'linear', slurmBot)
         #add_threshgrad_optimizations(cellName, 'exponential', slurmBot)
  
      print '%d jobs will be queued...' % len(slurmBot.queue)
      slurmBot.run_and_wait()

      exit()
   
   if cmd == 'missing':
   
      cellNamesStr = sys.argv[2]
      cellNames = cellNamesStr.split(',')

      rerun = False
      if len(sys.argv) > 3:
         rstr = sys.argv[3]
         rerun = rstr == 'run'
      
      missingJobBatchFiles = []
      for cellName in cellNames:
         jbfs = find_missing_job_outputs(cellName, 'strflab.stft.LARS.nlType_linear.lambda_%0.3f.JN_%d.txt', lambda2Vals, 3, 'lars.nlType_linear.lambda_%0.3f.JN_%d')
         missingJobBatchFiles.extend(jbfs)
         jbfs = find_missing_job_outputs(cellName, 'strflab.stft.LARS.nlType_exponential.lambda_%0.3f.JN_%d.txt', lambda2Vals, 3, 'lars.nlType_exponential.lambda_%0.3f.JN_%d')
         missingJobBatchFiles.extend(jbfs)
         jbfs = find_missing_job_outputs(cellName, 'strflab.stft.threshgrad.nlType_linear.thresh_%0.2f.JN_%d.txt', threshVals, 3, 'threshgrad.nlType_linear.thresh_%0.2f.JN_%d')
         missingJobBatchFiles.extend(jbfs)
         jbfs = find_missing_job_outputs(cellName, 'strflab.stft.threshgrad.nlType_exponential.thresh_%0.2f.JN_%d.txt', threshVals, 3, 'threshgrad.nlType_exponential.thresh_%0.2f.JN_%d')
         missingJobBatchFiles.extend(jbfs)

      if rerun:
         for jbf in missingJobBatchFiles:
            print 'Running script: %s' % jbf
            jid = slurmtools.slurm_sbatch_from_file(jbf)
            print '\tJob ID: %d' % jid
            time.sleep(0.1)
         
      exit()

   

   if cmd == 'missingfinal':

      cellNamesStr = sys.argv[2]
      cellNames = cellNamesStr.split(',')

      rerun = False
      if len(sys.argv) > 3:
         rstr = sys.argv[3]
         rerun = rstr == 'run'
      
      missingJobBatchFiles = []
      for cellName in cellNames:
         sbname = '%s.threshgrad.best.nlType_%%s.sbatch' % cellName
         jbfs = find_missing_final_job_outputs(cellName, 'strflab.stft.threshgrad.best.nlType_%s.txt', sbname)
         missingJobBatchFiles.extend(jbfs)
         sbname = '%s.lars.best.nlType_%%s.sbatch' % cellName
         jbfs = find_missing_final_job_outputs(cellName, 'strflab.stft.LARS.best.nlType_%s.txt', sbname)
         missingJobBatchFiles.extend(jbfs)
         sbname = '%s.directfit.nlType_%%s.sbatch' % cellName
         jbfs = find_missing_final_job_outputs(cellName, 'strflab.stft.directfit.nlType_%s.txt', sbname)
         missingJobBatchFiles.extend(jbfs)

      exit()

   
   if cmd == 'rename':

      """
      cellNamesStr = sys.argv[2]
      cellNames = cellNamesStr.split(',')
      
      for cellName in cellNames:
         print 'Renaming output for %s...' % cellName
         #rename_job_outputs(cellName, 'strflab.stft.lin.LARS.nlType_linear.lambda_%0.3f.JN_%d', 'strflab.stft.LARS.nlType_linear.lambda_%0.3f.JN_%d', lambda2Vals, 3)
         #rename_job_outputs(cellName, 'strflab.stft.lin.LARS.nlType_exponential.lambda_%0.3f.JN_%d', 'strflab.stft.LARS.nlType_exponential.lambda_%0.3f.JN_%d', lambda2Vals, 3)
         rename_job_outputs(cellName, 'strflab.stft.linear.threshgrad.nlType_linear.thresh_%0.2f.JN_%d', 'strflab.stft.threshgrad.nlType_linear.thresh_%0.2f.JN_%d', threshVals, 3)
         rename_job_outputs(cellName, 'strflab.stft.exponential.threshgrad.nlType_exponential.thresh_%0.2f.JN_%d', 'strflab.stft.threshgrad.nlType_exponential.thresh_%0.2f.JN_%d', threshVals, 3)
         
      exit()
      """
